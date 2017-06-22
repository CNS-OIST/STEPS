###___license_placeholder___###

from steps_wmrk4 cimport *
from steps_wmdirect cimport *
from steps_tetexact cimport *
from steps_tetode cimport *
from steps_solver cimport *

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

    def __init__(self, _py_Model m, _py_Geom g, _py_RNG r):
        self._ptr = new Wmrk4(m.ptr(), g.ptr(), r.ptr() )
        _py_API.__init__(self, m, g, r)

    def getSolverName(self, ):
        return self.ptrx().getSolverName()

    def getSolverDesc(self, ):
        return self.ptrx().getSolverDesc()

    def getSolverAuthors(self, ):
        return self.ptrx().getSolverAuthors()

    def getSolverEmail(self, ):
        return self.ptrx().getSolverEmail()

    def reset(self, ):
        self.ptrx().reset()

    def run(self, double endtime):
        self.ptrx().run(endtime)

    def advance(self, double adv):
        self.ptrx().advance(adv)

    def step(self, ):
        self.ptrx().step()

    def setDT(self, double dt):
        self.ptrx().setDT(dt)

    def setRk4DT(self, double dt):
        self.ptrx().setRk4DT(dt)

    def getTime(self, ):
        return self.ptrx().getTime()

    def checkpoint(self, std.string file_name):
        self.ptrx().checkpoint(file_name)

    def restore(self, std.string file_name):
        self.ptrx().restore(file_name)

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
        self._ptr = new Wmdirect( m.ptr(), g.ptr(), r.ptr() )
        #super(self.__class__, self).__init__(m,g,r)
        _py_API.__init__(self, m, g, r)

    def checkpoint(self, std.string file_name):
        self.ptrd().checkpoint(file_name)

    def restore(self, std.string file_name):
        self.ptrd().restore(file_name)

    def getSolverName(self, ):
        return self.ptrd().getSolverName()

    def getSolverDesc(self, ):
        return self.ptrd().getSolverDesc()

    def getSolverAuthors(self, ):
        return self.ptrd().getSolverAuthors()

    def getSolverEmail(self, ):
        return self.ptrd().getSolverEmail()

    def reset(self, ):
        self.ptrd().reset()

    def run(self, double endtime):
        self.ptrd().run(endtime)

    def advance(self, double adv):
        self.ptrd().advance(adv)

    def step(self, ):
        self.ptrd().step()

    def getTime(self, ):
        return self.ptrd().getTime()

    def getA0(self, ):
        return self.ptrd().getA0()

    def getNSteps(self, ):
        return self.ptr().getNSteps()

    def setTime(self, double time):
        self.ptrd().setTime(time)

    def setNSteps(self, unsigned int nsteps):
        self.ptrd().setNSteps(nsteps)

    # def addKProc(self, steps.wmdirect.KProc* kp):
    #     return _py_void.from_ref(self.ptr().addKProc(kp.ptr()))

    def countKProcs(self, ):
        return self.ptrd().countKProcs()

    @staticmethod
    cdef _py_Wmdirect from_ptr(Wmdirect *ptr):
        cdef _py_Wmdirect obj = _py_Wmdirect.__new__(_py_Wmdirect)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Wmdirect from_ref(const Wmdirect &ref):
        return _py_Wmdirect.from_ptr(<Wmdirect*>&ref)


# ======================================================================================================================
# Python bindings to namespace steps::tetexact
# ======================================================================================================================

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TEKProc(_py__base):
    "Python wrapper class for TEKProc (steps::tetexact::KProc)"
# ----------------------------------------------------------------------------------------------------------------------
    cdef TEKProc * ptr(self):
        return <TEKProc* > self._ptr

    @staticmethod
    cdef _py_TEKProc from_ptr(TEKProc *ptr):
        cdef _py_TEKProc obj = _py_TEKProc.__new__(_py_TEKProc)
        obj._ptr = ptr
        return obj


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Tetexact(_py_API):
    "Python wrapper class for Tetexact"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Tetexact] _autodealoc
    cdef Tetexact *ptrx(self):
        return <Tetexact*> self._ptr

    def __init__(self, _py_Model m, _py_Geom g, _py_RNG r, int calcMembPot=0):
        self._ptr = new Tetexact(m.ptr(), g.ptr(), r.ptr(), calcMembPot)
        _py_API.__init__(self, m, g, r)

    def getSolverName(self, ):
        return self.ptrx().getSolverName()

    def getSolverDesc(self, ):
        return self.ptrx().getSolverDesc()

    def getSolverAuthors(self, ):
        return self.ptrx().getSolverAuthors()

    def getSolverEmail(self, ):
        return self.ptrx().getSolverEmail()

    def reset(self, ):
        self.ptrx().reset()

    def run(self, double endtime):
        self.ptrx().run(endtime)

    def advance(self, double adv):
        self.ptrx().advance(adv)

    def step(self, ):
        self.ptrx().step()

    def checkpoint(self, std.string file_name):
        self.ptrx().checkpoint(file_name)

    def restore(self, std.string file_name):
        self.ptrx().restore(file_name)

    def setEfieldDT(self, double efdt):
        self.ptrx().setEfieldDT(efdt)

    def efdt(self, ):
        return self.ptrx().efdt()

    def setTemp(self, double t):
        self.ptrx().setTemp(t)

    def getTemp(self, ):
        return self.ptrx().getTemp()

    def saveMembOpt(self, std.string opt_file_name):
        self.ptrx().saveMembOpt(opt_file_name)

    def getTime(self, ):
        return self.ptrx().getTime()

    def getA0(self, ):
        return self.ptrx().getA0()

    def getNSteps(self, ):
        return self.ptrx().getNSteps()

    def setTime(self, double time):
        self.ptrx().setTime(time)

    def setNSteps(self, unsigned int nsteps):
        self.ptrx().setNSteps(nsteps)

    def addKProc(self, _py_TEKProc kp):
        self.ptrx().addKProc(kp.ptr())

    def countKProcs(self, ):
        return self.ptrx().countKProcs()

    def mesh(self, ):
        return _py_Tetmesh.from_ptr(self.ptrx().mesh())

    @staticmethod
    cdef _py_Tetexact from_ptr(Tetexact *ptr):
        cdef _py_Tetexact obj = _py_Tetexact.__new__(_py_Tetexact)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Tetexact from_ref(const Tetexact &ref):
        return _py_Tetexact.from_ptr(<Tetexact*>&ref)


# ======================================================================================================================
# Python bindings to namespace steps::tetode
# ======================================================================================================================

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TetODE(_py_API):
    "Python wrapper class for TetODE"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[TetODE] _autodealoc
    cdef TetODE *ptrx(self):
        return <TetODE*> self._ptr

    def __init__(self, _py_Model m, _py_Geom g, _py_RNG r=None, int calcMembPot=0):
        self._ptr = new TetODE(m.ptr(), g.ptr(), r.ptr() if r else NULL, calcMembPot)
        _py_API.__init__(self, m, g, r)

    def getSolverName(self, ):
        return self.ptrx().getSolverName()

    def getSolverDesc(self, ):
        return self.ptrx().getSolverDesc()

    def getSolverAuthors(self, ):
        return self.ptrx().getSolverAuthors()

    def getSolverEmail(self, ):
        return self.ptrx().getSolverEmail()

    def checkpoint(self, std.string file_name):
        self.ptrx().checkpoint(file_name)

    def restore(self, std.string file_name):
        self.ptrx().restore(file_name)

    def getTime(self, ):
        return self.ptrx().getTime()

    def getTemp(self, ):
        return self.ptrx().getTemp()

    def setTemp(self, double t):
        return self.ptrx().setTemp(t)

    def reset(self, ):
        self.ptrx().reset()

    def run(self, double endtime):
        self.ptrx().run(endtime)

    def advance(self, double adv):
        self.ptrx().advance(adv)

    def mesh(self, ):
        return _py_Tetmesh.from_ptr(self.ptrx().mesh())

    def setTolerances(self, double atol, double rtol):
        self.ptrx().setTolerances(atol, rtol)

    def setMaxNumSteps(self, unsigned int maxn):
        self.ptrx().setMaxNumSteps(maxn)

    def efflag(self, ):
        return self.ptrx().efflag()

    def neftets(self, ):
        return self.ptrx().neftets()

    def neftris(self, ):
        return self.ptrx().neftris()

    def nefverts(self, ):
        return self.ptrx().nefverts()

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
    model = None
    geom = None

    #Constants
    EF_NONE      = steps_solver.EF_NONE
    EF_DEFAULT   = steps_solver.EF_DEFAULT
    EF_DV_BDSYS  = steps_solver.EF_DV_BDSYS
    EF_DV_SLUSYS = steps_solver.EF_DV_SLUSYS
    EF_DV_PETSC  = steps_solver.EF_DV_PETSC

    cdef API *ptr(self):
        return <API*> self._ptr

    # ---- VIRTUAL - doesnt call original constructor ------
    def __init__(self, _py_Model m, _py_Geom g, _py_RNG r):
        self.model = m
        self.geom = g

    def getSolverName(self, ):
        return self.ptr().getSolverName()

    def getSolverDesc(self, ):
        return self.ptr().getSolverDesc()

    def getSolverAuthors(self, ):
        return self.ptr().getSolverAuthors()

    def getSolverEmail(self, ):
        return self.ptr().getSolverEmail()

    def checkpoint(self, std.string file_name):
        self.ptr().checkpoint(file_name)

    def restore(self, std.string file_name):
        self.ptr().restore(file_name)

    def reset(self, ):
        self.ptr().reset()

    def run(self, double endtime):
        self.ptr().run(endtime)

    def advance(self, double adv):
        self.ptr().advance(adv)

    def step(self, ):
        self.ptr().step()

    def setRk4DT(self, double dt):
        self.ptr().setRk4DT(dt)

    def setDT(self, double dt):
        self.ptr().setDT(dt)

    def setEfieldDT(self, double efdt):
        self.ptr().setEfieldDT(efdt)

    def setNSteps(self, unsigned int nsteps):
        self.ptr().setNSteps(nsteps)

    def setTime(self, double time):
        self.ptr().setTime(time)

    def setTemp(self, double temp):
        self.ptr().setTemp(temp)

    def getTime(self, ):
        return self.ptr().getTime()

    def getRk4DT(self, ):
        return self.ptr().getRk4DT()

    def getDT(self, ):
        return self.ptr().getDT()

    def getEfieldDT(self, ):
        return self.ptr().getEfieldDT()

    def getTemp(self, ):
        return self.ptr().getTemp()

    def getA0(self, ):
        return self.ptr().getA0()

    def getNSteps(self, ):
        return self.ptr().getNSteps()

    def getCompVol(self, std.string c):
        return self.ptr().getCompVol(c)

    def setCompVol(self, std.string c, double vol):
        self.ptr().setCompVol(c, vol)

    def getCompCount(self, std.string c, std.string s):
        return self.ptr().getCompCount(c, s)

    def setCompCount(self, std.string c, std.string s, double n):
        self.ptr().setCompCount(c, s, n)

    def getCompAmount(self, std.string c, std.string s):
        return self.ptr().getCompAmount(c, s)

    def setCompAmount(self, std.string c, std.string s, double a):
        self.ptr().setCompAmount(c, s, a)

    def getCompConc(self, std.string c, std.string s):
        return self.ptr().getCompConc(c, s)

    def setCompConc(self, std.string c, std.string s, double conc):
        self.ptr().setCompConc(c, s, conc)

    def getCompClamped(self, std.string c, std.string s):
        return self.ptr().getCompClamped(c, s)

    def setCompClamped(self, std.string c, std.string s, bool b):
        self.ptr().setCompClamped(c, s, b)

    def getCompReacK(self, std.string c, std.string r):
        return self.ptr().getCompReacK(c, r)

    def setCompReacK(self, std.string c, std.string r, double kf):
        self.ptr().setCompReacK(c, r, kf)

    def getCompReacActive(self, std.string c, std.string r):
        return self.ptr().getCompReacActive(c, r)

    def setCompReacActive(self, std.string c, std.string r, bool a):
        self.ptr().setCompReacActive(c, r, a)

    def getCompDiffD(self, std.string c, std.string d):
        return self.ptr().getCompDiffD(c, d)

    def setCompDiffD(self, std.string c, std.string d, double dcst):
        self.ptr().setCompDiffD(c, d, dcst)

    def getCompDiffActive(self, std.string c, std.string d):
        return self.ptr().getCompDiffActive(c, d)

    def setCompDiffActive(self, std.string c, std.string d, bool act):
        self.ptr().setCompDiffActive(c, d, act)

    def getCompReacC(self, std.string c, std.string r):
        return self.ptr().getCompReacC(c, r)

    def getCompReacH(self, std.string c, std.string r):
        return self.ptr().getCompReacH(c, r)

    def getCompReacA(self, std.string c, std.string r):
        return self.ptr().getCompReacA(c, r)

    def getCompReacExtent(self, std.string c, std.string r):
        return self.ptr().getCompReacExtent(c, r)

    def resetCompReacExtent(self, std.string c, std.string r):
        self.ptr().resetCompReacExtent(c, r)

    def getTetVol(self, unsigned int tidx):
        return self.ptr().getTetVol(tidx)

    def setTetVol(self, unsigned int tidx, double vol):
        self.ptr().setTetVol(tidx, vol)

    def getTetSpecDefined(self, unsigned int tidx, std.string s):
        return self.ptr().getTetSpecDefined(tidx, s)

    def getTetCount(self, unsigned int tidx, std.string s):
        return self.ptr().getTetCount(tidx, s)

    def setTetCount(self, unsigned int tidx, std.string s, double n):
        self.ptr().setTetCount(tidx, s, n)

    def getTetAmount(self, unsigned int tidx, std.string s):
        return self.ptr().getTetAmount(tidx, s)

    def setTetAmount(self, unsigned int tidx, std.string s, double m):
        self.ptr().setTetAmount(tidx, s, m)

    def getTetConc(self, unsigned int tidx, std.string s):
        return self.ptr().getTetConc(tidx, s)

    def setTetConc(self, unsigned int tidx, std.string s, double c):
        self.ptr().setTetConc(tidx, s, c)

    def getTetClamped(self, unsigned int tidx, std.string s):
        return self.ptr().getTetClamped(tidx, s)

    def setTetClamped(self, unsigned int tidx, std.string s, bool buf):
        self.ptr().setTetClamped(tidx, s, buf)

    def getTetReacK(self, unsigned int tidx, std.string r):
        return self.ptr().getTetReacK(tidx, r)

    def setTetReacK(self, unsigned int tidx, std.string r, double kf):
        self.ptr().setTetReacK(tidx, r, kf)

    def getTetReacActive(self, unsigned int tidx, std.string r):
        return self.ptr().getTetReacActive(tidx, r)

    def setTetReacActive(self, unsigned int tidx, std.string r, bool act):
        self.ptr().setTetReacActive(tidx, r, act)

    def getTetDiffD(self, unsigned int tidx, std.string d, unsigned int direction_tet=numeric_limits[uint].max()):
        return self.ptr().getTetDiffD(tidx, d, direction_tet)

    def setTetDiffD(self, unsigned int tidx, std.string d, double dk, unsigned int direction_tet=numeric_limits[uint].max()):
        self.ptr().setTetDiffD(tidx, d, dk, direction_tet)

    def getTetDiffActive(self, unsigned int tidx, std.string d):
        return self.ptr().getTetDiffActive(tidx, d)

    def setTetDiffActive(self, unsigned int tidx, std.string d, bool act):
        self.ptr().setTetDiffActive(tidx, d, act)

    def getTetReacC(self, unsigned int tidx, std.string r):
        return self.ptr().getTetReacC(tidx, r)

    def getTetReacH(self, unsigned int tidx, std.string r):
        return self.ptr().getTetReacH(tidx, r)

    def getTetReacA(self, unsigned int tidx, std.string r):
        return self.ptr().getTetReacA(tidx, r)

    def getTetDiffA(self, unsigned int tidx, std.string d):
        return self.ptr().getTetDiffA(tidx, d)

    def getTetV(self, unsigned int tidx):
        return self.ptr().getTetV(tidx)

    def setTetV(self, unsigned int tidx, double v):
        self.ptr().setTetV(tidx, v)

    def getTetVClamped(self, unsigned int tidx):
        return self.ptr().getTetVClamped(tidx)

    def setTetVClamped(self, unsigned int tidx, bool cl):
        self.ptr().setTetVClamped(tidx, cl)

    def getPatchArea(self, std.string p):
        return self.ptr().getPatchArea(p)

    def setPatchArea(self, std.string p, double area):
        self.ptr().setPatchArea(p, area)

    def getPatchCount(self, std.string p, std.string s):
        return self.ptr().getPatchCount(p, s)

    def setPatchCount(self, std.string p, std.string s, double n):
        self.ptr().setPatchCount(p, s, n)

    def getPatchAmount(self, std.string p, std.string s):
        return self.ptr().getPatchAmount(p, s)

    def setPatchAmount(self, std.string p, std.string s, double a):
        self.ptr().setPatchAmount(p, s, a)

    def getPatchClamped(self, std.string p, std.string s):
        return self.ptr().getPatchClamped(p, s)

    def setPatchClamped(self, std.string p, std.string s, bool buf):
        self.ptr().setPatchClamped(p, s, buf)

    def getPatchSReacK(self, std.string p, std.string r):
        return self.ptr().getPatchSReacK(p, r)

    def setPatchSReacK(self, std.string p, std.string r, double kf):
        self.ptr().setPatchSReacK(p, r, kf)

    def getPatchSReacActive(self, std.string p, std.string r):
        return self.ptr().getPatchSReacActive(p, r)

    def setPatchSReacActive(self, std.string p, std.string r, bool a):
        self.ptr().setPatchSReacActive(p, r, a)

    def getPatchSReacC(self, std.string p, std.string r):
        return self.ptr().getPatchSReacC(p, r)

    def getPatchSReacH(self, std.string p, std.string r):
        return self.ptr().getPatchSReacH(p, r)

    def getPatchSReacA(self, std.string p, std.string r):
        return self.ptr().getPatchSReacA(p, r)

    def getPatchSReacExtent(self, std.string p, std.string r):
        return self.ptr().getPatchSReacExtent(p, r)

    def resetPatchSReacExtent(self, std.string p, std.string r):
        self.ptr().resetPatchSReacExtent(p, r)

    def getPatchVDepSReacActive(self, std.string p, std.string vsr):
        return self.ptr().getPatchVDepSReacActive(p, vsr)

    def setPatchVDepSReacActive(self, std.string p, std.string vsr, bool a):
        self.ptr().setPatchVDepSReacActive(p, vsr, a)

    def setDiffBoundaryDiffusionActive(self, std.string db, std.string s, bool act):
        self.ptr().setDiffBoundaryDiffusionActive(db, s, act)

    def getDiffBoundaryDiffusionActive(self, std.string db, std.string s):
        return self.ptr().getDiffBoundaryDiffusionActive(db, s)

    def setDiffBoundaryDcst(self, std.string db, std.string s, double dcst, std.string direction_comp=""):
        self.ptr().setDiffBoundaryDcst(db, s, dcst, direction_comp)

    def setSDiffBoundaryDiffusionActive(self, std.string db, std.string s, bool act):
        self.ptr().setSDiffBoundaryDiffusionActive(db, s, act)

    def getSDiffBoundaryDiffusionActive(self, std.string db, std.string s):
        return self.ptr().getSDiffBoundaryDiffusionActive(db, s)

    def setSDiffBoundaryDcst(self, std.string db, std.string s, double dcst, std.string direction_patch=""):
        self.ptr().setSDiffBoundaryDcst(db, s, dcst, direction_patch)

    def getTriArea(self, unsigned int tidx):
        return self.ptr().getTriArea(tidx)

    def setTriArea(self, unsigned int tidx, double area):
        self.ptr().setTriArea(tidx, area)

    def getTriSpecDefined(self, unsigned int tidx, std.string s):
        return self.ptr().getTriSpecDefined(tidx, s)

    def getTriCount(self, unsigned int tidx, std.string s):
        return self.ptr().getTriCount(tidx, s)

    def setTriCount(self, unsigned int tidx, std.string s, double n):
        self.ptr().setTriCount(tidx, s, n)

    def getTriAmount(self, unsigned int tidx, std.string s):
        return self.ptr().getTriAmount(tidx, s)

    def setTriAmount(self, unsigned int tidx, std.string s, double m):
        self.ptr().setTriAmount(tidx, s, m)

    def getTriClamped(self, unsigned int tidx, std.string s):
        return self.ptr().getTriClamped(tidx, s)

    def setTriClamped(self, unsigned int tidx, std.string s, bool buf):
        self.ptr().setTriClamped(tidx, s, buf)

    def getTriSReacK(self, unsigned int tidx, std.string r):
        return self.ptr().getTriSReacK(tidx, r)

    def setTriSReacK(self, unsigned int tidx, std.string r, double kf):
        self.ptr().setTriSReacK(tidx, r, kf)

    def getTriSReacActive(self, unsigned int tidx, std.string r):
        return self.ptr().getTriSReacActive(tidx, r)

    def setTriSReacActive(self, unsigned int tidx, std.string r, bool act):
        self.ptr().setTriSReacActive(tidx, r, act)

    def getTriSReacC(self, unsigned int tidx, std.string r):
        return self.ptr().getTriSReacC(tidx, r)

    def getTriSReacH(self, unsigned int tidx, std.string r):
        return self.ptr().getTriSReacH(tidx, r)

    def getTriSReacA(self, unsigned int tidx, std.string r):
        return self.ptr().getTriSReacA(tidx, r)

    def getTriDiffD(self, unsigned int tidx, std.string d, unsigned int direction_tri=numeric_limits[uint].max()):
        return self.ptr().getTriDiffD(tidx, d, direction_tri)

    def getTriSDiffD(self, unsigned int tidx, std.string d, unsigned int direction_tri=numeric_limits[uint].max()):
        return self.ptr().getTriSDiffD(tidx, d, direction_tri)

    def setTriDiffD(self, unsigned int tidx, std.string d, double dk, unsigned int direction_tri=numeric_limits[uint].max()):
        self.ptr().setTriDiffD(tidx, d, dk, direction_tri)

    def setTriSDiffD(self, unsigned int tidx, std.string d, double dk, unsigned int direction_tri=numeric_limits[uint].max()):
        self.ptr().setTriSDiffD(tidx, d, dk, direction_tri)

    def getTriV(self, unsigned int tidx):
        return self.ptr().getTriV(tidx)

    def setTriV(self, unsigned int tidx, double v):
        self.ptr().setTriV(tidx, v)

    def getTriVClamped(self, unsigned int tidx):
        return self.ptr().getTriVClamped(tidx)

    def setTriVClamped(self, unsigned int tidx, bool cl):
        self.ptr().setTriVClamped(tidx, cl)

    def getTriOhmicI(self, unsigned int tidx, std.string oc=''):
        if oc == '':
            return self.ptr().getTriOhmicI(tidx)
        return self.ptr().getTriOhmicI(tidx, oc)

    def getTriGHKI(self, unsigned int tidx, std.string ghk=''):
        if ghk == '':
            return self.ptr().getTriGHKI(tidx)
        return self.ptr().getTriGHKI(tidx, ghk)

    def getTriI(self, unsigned int tidx):
        return self.ptr().getTriI(tidx)

    def setTriIClamp(self, unsigned int tidx, double i):
        self.ptr().setTriIClamp(tidx, i)

    def getTriVDepSReacActive(self, unsigned int tidx, std.string vsr):
        return self.ptr().getTriVDepSReacActive(tidx, vsr)

    def setTriVDepSReacActive(self, unsigned int tidx, std.string vsr, bool act):
        self.ptr().setTriVDepSReacActive(tidx, vsr, act)

    def setTriCapac(self, unsigned int tidx, double cm):
        self.ptr().setTriCapac(tidx, cm)

    def getVertV(self, unsigned int vidx):
        return self.ptr().getVertV(vidx)

    def setVertV(self, unsigned int vidx, double v):
        self.ptr().setVertV(vidx, v)

    def getVertVClamped(self, unsigned int vidx):
        return self.ptr().getVertVClamped(vidx)

    def setVertVClamped(self, unsigned int vidx, bool cl):
        self.ptr().setVertVClamped(vidx, cl)

    def setVertIClamp(self, unsigned int vidx, double i):
        self.ptr().setVertIClamp(vidx, i)

    def setMembPotential(self, std.string m, double v):
        self.ptr().setMembPotential(m, v)

    def setMembCapac(self, std.string m, double cm):
        self.ptr().setMembCapac(m, cm)

    def setMembVolRes(self, std.string m, double ro):
        self.ptr().setMembVolRes(m, ro)

    def setMembRes(self, std.string m, double ro, double vrev):
        self.ptr().setMembRes(m, ro, vrev)

    def getNComps(self, ):
        return self.ptr().getNComps()

    def getNPatches(self, ):
        return self.ptr().getNPatches()

    def getCompName(self, unsigned int c_idx):
        return self.ptr().getCompName(c_idx)

    def getPatchName(self, unsigned int p_idx):
        return self.ptr().getPatchName(p_idx)

    def getNCompSpecs(self, unsigned int c_idx):
        return self.ptr().getNCompSpecs(c_idx)

    def getNPatchSpecs(self, unsigned int p_idx):
        return self.ptr().getNPatchSpecs(p_idx)

    def getCompSpecName(self, unsigned int c_idx, unsigned int s_idx):
        return self.ptr().getCompSpecName(c_idx, s_idx)

    def getPatchSpecName(self, unsigned int p_idx, unsigned int s_idx):
        return self.ptr().getPatchSpecName(p_idx, s_idx)

    def getBatchTetCounts(self, std.vector[unsigned int] tets, std.string s):
        return self.ptr().getBatchTetCounts(tets, s)

    def getBatchTriCounts(self, std.vector[unsigned int] tris, std.string s):
        return self.ptr().getBatchTriCounts(tris, s)

    def getBatchTetCountsNP(self, uint[:] indices, std.string s, double[:] counts):
        self.ptr().getBatchTetCountsNP(&indices[0], indices.shape[0], s, &counts[0], counts.shape[0])

    def getBatchTriCountsNP(self, uint[:] indices, std.string s, double[:] counts):
        self.ptr().getBatchTriCountsNP(&indices[0], indices.shape[0], s, &counts[0], counts.shape[0])

    def getROITetCounts(self, std.string ROI_id, std.string s):
        return self.ptr().getROITetCounts(ROI_id, s)

    def getROITriCounts(self, std.string ROI_id, std.string s):
        return self.ptr().getROITriCounts(ROI_id, s)

    def getROITetCountsNP(self, std.string ROI_id, std.string s, double[:] counts):
        self.ptr().getROITetCountsNP(ROI_id, s, &counts[0], counts.shape[0])

    def getROITriCountsNP(self, std.string ROI_id, std.string s, double[:] counts):
        self.ptr().getROITriCountsNP(ROI_id, s, &counts[0], counts.shape[0])

    def getROIVol(self, std.string ROI_id):
        return self.ptr().getROIVol(ROI_id)

    def getROIArea(self, std.string ROI_id):
        return self.ptr().getROIArea(ROI_id)

    def getROICount(self, std.string ROI_id, std.string s):
        return self.ptr().getROICount(ROI_id, s)

    def setROICount(self, std.string ROI_id, std.string s, double count):
        self.ptr().setROICount(ROI_id, s, count)

    def getROIAmount(self, std.string ROI_id, std.string s):
        return self.ptr().getROIAmount(ROI_id, s)

    def getROIConc(self, std.string ROI_id, std.string s):
        return self.ptr().getROIConc(ROI_id, s)

    def setROIConc(self, std.string ROI_id, std.string s, double conc):
        self.ptr().setROIConc(ROI_id, s, conc)

    def setROIClamped(self, std.string ROI_id, std.string s, bool b):
        self.ptr().setROIClamped(ROI_id, s, b)

    def setROIReacK(self, std.string ROI_id, std.string r, double kf):
        self.ptr().setROIReacK(ROI_id, r, kf)

    def setROISReacK(self, std.string ROI_id, std.string sr, double kf):
        self.ptr().setROISReacK(ROI_id, sr, kf)

    def setROIDiffD(self, std.string ROI_id, std.string d, double dk):
        self.ptr().setROIDiffD(ROI_id, d, dk)

    def setROIReacActive(self, std.string ROI_id, std.string r, bool a):
        self.ptr().setROIReacActive(ROI_id, r, a)

    def setROISReacActive(self, std.string ROI_id, std.string sr, bool a):
        self.ptr().setROISReacActive(ROI_id, sr, a)

    def setROIDiffActive(self, std.string ROI_id, std.string d, bool act):
        self.ptr().setROIDiffActive(ROI_id, d, act)

    def setROIVDepSReacActive(self, std.string ROI_id, std.string vsr, bool a):
        self.ptr().setROIVDepSReacActive(ROI_id, vsr, a)

    def getROIReacExtent(self, std.string ROI_id, std.string r):
        return self.ptr().getROIReacExtent(ROI_id, r)

    def resetROIReacExtent(self, std.string ROI_id, std.string r):
        self.ptr().resetROIReacExtent(ROI_id, r)

    def getROISReacExtent(self, std.string ROI_id, std.string sr):
        return self.ptr().getROISReacExtent(ROI_id, sr)

    def resetROISReacExtent(self, std.string ROI_id, std.string sr):
        self.ptr().resetROISReacExtent(ROI_id, sr)

    def getROIDiffExtent(self, std.string ROI_id, std.string d):
        return self.ptr().getROIDiffExtent(ROI_id, d)

    def resetROIDiffExtent(self, std.string ROI_id, std.string s):
        self.ptr().resetROIDiffExtent(ROI_id, s)

    @staticmethod
    cdef _py_API from_ptr(API *ptr):
        cdef _py_API obj = _py_API.__new__(_py_API )
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_API from_ref(const API &ref):
        return _py_API.from_ptr(<API*>&ref)
