###___license_placeholder___###

include "cysteps.pyx"

# ======================================================================================================================
# Python bindings to namespace steps::mpi
# ======================================================================================================================
cimport steps_mpi
from steps_mpi cimport TetOpSplitP

def mpiInit():
    steps_mpi.mpiInit()
    
def getRank():
    return steps_mpi.getRank()

def getNHosts():
    return steps_mpi.getNHosts()

def mpiFinish():
    steps_mpi.mpiFinish()
    
# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TetOpSplitP(_py_API):
    """Bindings for MPI TetOpSplitP"""
# ----------------------------------------------------------------------------------------------------------------------
    cdef TetOpSplitP *ptrx(self):
        return <TetOpSplitP*> self._ptr
    
    def __init__(self, _py_Model model, _py_Geom geom, _py_RNG rng, int calcMembPot=0, std.vector[uint] tet_hosts = [], dict tri_hosts = {}, std.vector[uint] wm_hosts = []):
        cdef std.map[uint, uint] _tri_hosts
        for key, elem in tri_hosts.items():
            _tri_hosts[key] = elem
        # We constructed a map. Now call constructor
        self._ptr = new TetOpSplitP(model.ptr(), geom.ptr(), rng.ptr(), calcMembPot, tet_hosts, _tri_hosts, wm_hosts)
        
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

    ### All these entities are specific to steps::mpi. We would need to wrapp them all ###
    # def addKProc(self, steps.mpi.tetopsplit.KProc* kp, int host):
    #     return PY_unsigned int.from_ref(self.ptr().addKProc(kp.ptr(), deref(host.ptr())))
    #
    # def addDiff(self, steps.mpi.tetopsplit.Diff* diff):
    #     return PY_void.from_ref(self.ptr().addDiff(diff.ptr()))
    #
    # def addSDiff(self, steps.mpi.tetopsplit.SDiff* sdiff):
    #     return PY_void.from_ref(self.ptr().addSDiff(sdiff.ptr()))
    #
    # def patches(self, ):
    #     return PY_std.vector[Patch*].from_ref(self.ptr().patches())

    def countKProcs(self, ):
        return self.ptrx().countKProcs()

    def mesh(self, ):
        return _py_Tetmesh.from_ptr(self.ptrx().mesh())

    def a0(self, ):
        return self.ptxr().a0()

    def efflag(self, ):
        return self.ptrx().efflag()

    def neftets(self, ):
        return self.ptrx().neftets()

    def neftris(self, ):
        return self.ptrx().neftris()

    def nefverts(self, ):
        return self.ptrx().nefverts()

    def getBatchTetCounts(self, std.vector[uint] tets, std.string s):
        return self.ptrx().getBatchTetCounts(tets, s)

    def getBatchTriCounts(self, std.vector[uint] tris, std.string s):
        return self.ptrx().getBatchTriCounts(tris, s)

    # ---------------------------------------------------------------------------------
    # NUMPY section - we accept numpy arrays and generically typed memory-views
    # ---------------------------------------------------------------------------------
    def getBatchTetCountsNP(self, uint[:] index_array, std.string s, double[:] counts):
        self.ptrx().getBatchTetCountsNP(&index_array[0], index_array.shape[0], s, &counts[0], counts.shape[0])

    def getBatchTriCountsNP(self, uint[:] index_array, std.string s, double[:] counts):
        self.ptrx().getBatchTriCountsNP(&index_array[0], index_array.shape[0], s, &counts[0], counts.shape[0])

    def sumBatchTetCountsNP(self, uint[:] tet_array, std.string s):
        return self.ptrx().sumBatchTetCountsNP(&tet_array[0], tet_array.shape[0], s)

    def sumBatchTriCountsNP(self, uint[:] tri_array, std.string s):
        return self.ptrx().sumBatchTriCountsNP(&tri_array[0], tri_array.shape[0], s)

    def sumBatchTriGHKIsNP(self, uint[:] tri_array, std.string ghk):
        return self.ptrx().sumBatchTriGHKIsNP(&tri_array[0], tri_array.shape[0], ghk)

    def sumBatchTriOhmicIsNP(self, uint[:] tri_array, std.string ghk):
        return self.ptrx().sumBatchTriOhmicIsNP(&tri_array[0], tri_array.shape[0], ghk)

    # ---------------------------------------------------------------------------------
    # ROI section
    # ---------------------------------------------------------------------------------
    # def getROITetCounts(self, std.string ROI_id, std.string s):
    #     return PY_std.vector[double].from_ref(self.ptr().getROITetCounts(deref(ROI_id.ptr()), deref(s.ptr())))
    #
    # def getROITriCounts(self, std.string ROI_id, std.string s):
    #     return PY_std.vector[double].from_ref(self.ptr().getROITriCounts(deref(ROI_id.ptr()), deref(s.ptr())))
    #
    # def getROITetCountsNP(self, std.string ROI_id, std.string s, double* counts, int output_size):
    #     return PY_void.from_ref(self.ptr().getROITetCountsNP(deref(ROI_id.ptr()), deref(s.ptr()), counts.ptr(), deref(output_size.ptr())))
    #
    # def getROITriCountsNP(self, std.string ROI_id, std.string s, double* counts, int output_size):
    #     return PY_void.from_ref(self.ptr().getROITriCountsNP(deref(ROI_id.ptr()), deref(s.ptr()), counts.ptr(), deref(output_size.ptr())))
    #
    # def getROIVol(self, std.string ROI_id):
    #     return PY_double.from_ref(self.ptr().getROIVol(deref(ROI_id.ptr())))
    #
    def getROIArea(self, std.string ROI_id):
        return self.ptr().getROIArea(ROI_id)
    #
    # def getROICount(self, std.string ROI_id, std.string s):
    #     return PY_double.from_ref(self.ptr().getROICount(deref(ROI_id.ptr()), deref(s.ptr())))
    #
    # def setROICount(self, std.string ROI_id, std.string s, double count):
    #     return PY_void.from_ref(self.ptr().setROICount(deref(ROI_id.ptr()), deref(s.ptr()), deref(count.ptr())))
    #
    # def getROIAmount(self, std.string ROI_id, std.string s):
    #     return PY_double.from_ref(self.ptr().getROIAmount(deref(ROI_id.ptr()), deref(s.ptr())))
    #
    # def getROIConc(self, std.string ROI_id, std.string s):
    #     return PY_double.from_ref(self.ptr().getROIConc(deref(ROI_id.ptr()), deref(s.ptr())))
    #
    # def setROIConc(self, std.string ROI_id, std.string s, double conc):
    #     return PY_void.from_ref(self.ptr().setROIConc(deref(ROI_id.ptr()), deref(s.ptr()), deref(conc.ptr())))
    #
    # def setROIClamped(self, std.string ROI_id, std.string s, bool b):
    #     return PY_void.from_ref(self.ptr().setROIClamped(deref(ROI_id.ptr()), deref(s.ptr()), deref(b.ptr())))
    #
    # def setROIReacK(self, std.string ROI_id, std.string r, double kf):
    #     return PY_void.from_ref(self.ptr().setROIReacK(deref(ROI_id.ptr()), deref(r.ptr()), deref(kf.ptr())))
    #
    # def setROISReacK(self, std.string ROI_id, std.string sr, double kf):
    #     return PY_void.from_ref(self.ptr().setROISReacK(deref(ROI_id.ptr()), deref(sr.ptr()), deref(kf.ptr())))
    #
    # def setROIDiffD(self, std.string ROI_id, std.string d, double dk):
    #     return PY_void.from_ref(self.ptr().setROIDiffD(deref(ROI_id.ptr()), deref(d.ptr()), deref(dk.ptr())))
    #
    # def setROIReacActive(self, std.string ROI_id, std.string r, bool a):
    #     return PY_void.from_ref(self.ptr().setROIReacActive(deref(ROI_id.ptr()), deref(r.ptr()), deref(a.ptr())))
    #
    # def setROISReacActive(self, std.string ROI_id, std.string sr, bool a):
    #     return PY_void.from_ref(self.ptr().setROISReacActive(deref(ROI_id.ptr()), deref(sr.ptr()), deref(a.ptr())))
    #
    # def setROIDiffActive(self, std.string ROI_id, std.string d, bool a):
    #     return PY_void.from_ref(self.ptr().setROIDiffActive(deref(ROI_id.ptr()), deref(d.ptr()), deref(a.ptr())))
    #
    # def setROIVDepSReacActive(self, std.string ROI_id, std.string vsr, bool a):
    #     return PY_void.from_ref(self.ptr().setROIVDepSReacActive(deref(ROI_id.ptr()), deref(vsr.ptr()), deref(a.ptr())))
    #
    # def getROIReacExtent(self, std.string ROI_id, std.string r):
    #     return PY_unsigned int.from_ref(self.ptr().getROIReacExtent(deref(ROI_id.ptr()), deref(r.ptr())))
    #
    # def resetROIReacExtent(self, std.string ROI_id, std.string r):
    #     return PY_void.from_ref(self.ptr().resetROIReacExtent(deref(ROI_id.ptr()), deref(r.ptr())))
    #
    # def getROISReacExtent(self, std.string ROI_id, std.string sr):
    #     return PY_unsigned int.from_ref(self.ptr().getROISReacExtent(deref(ROI_id.ptr()), deref(sr.ptr())))
    #
    # def resetROISReacExtent(self, std.string ROI_id, std.string sr):
    #     return PY_void.from_ref(self.ptr().resetROISReacExtent(deref(ROI_id.ptr()), deref(sr.ptr())))
    #
    # def getROIDiffExtent(self, std.string ROI_id, std.string d):
    #     return PY_unsigned int.from_ref(self.ptr().getROIDiffExtent(deref(ROI_id.ptr()), deref(d.ptr())))
    #
    # def resetROIDiffExtent(self, std.string ROI_id, std.string d):
    #     return PY_void.from_ref(self.ptr().resetROIDiffExtent(deref(ROI_id.ptr()), deref(d.ptr())))
    # ------------------------------------------------------------------------------------------------------------

    def setDiffApplyThreshold(self, int threshold):
        self.ptrx().setDiffApplyThreshold(threshold)

    def getTetHostRank(self, unsigned int tidx):
        return self.ptrx().getTetHostRank(tidx)

    def getTriHostRank(self, unsigned int tidx):
        return self.ptrx().getTriHostRank(tidx)

    def getWMVolHostRank(self, unsigned int idx):
        return self.ptrx().getWMVolHostRank(idx)

    def addNeighHost(self, int host):
        self.ptrx().addNeighHost(host)

#    def registerRemoteChange(self, steps_mpi.SubVolType svol_type, unsigned int idx, unsigned int slidx, int change):
#        self.ptrx().registerRemoteChange(svol_type, idx, slidx, change)

    def getReacExtent(self, bool local=False):
        return self.ptrx().getReacExtent(local)

    def getDiffExtent(self, bool local=False):
        return self.ptrx().getDiffExtent(local)

    def getNIteration(self, ):
        return self.ptrx().getNIteration()

#    def getNInternalDiffs(self, bool local=False):
#        return self.ptrx().getNInternalDiffs(local)

#    def getNExternalDiffs(self, bool local=False):
#        return self.ptrx().getNExternalDiffs(local)

    def getUpdPeriod(self, ):
        return self.ptrx().getUpdPeriod()

    def getCompTime(self, ):
        return self.ptrx().getCompTime()

    def getSyncTime(self, ):
        return self.ptrx().getSyncTime()

    def getIdleTime(self, ):
        return self.ptrx().getIdleTime()

    def getEFieldTime(self, ):
        return self.ptrx().getEFieldTime()

    def getRDTime(self, ):
        return self.ptrx().getRDTime()

    def getDataExchangeTime(self, ):
        return self.ptrx().getDataExchangeTime()

    @staticmethod
    cdef _py_TetOpSplitP from_ptr(TetOpSplitP *ptr):
        cdef _py_TetOpSplitP obj = _py_TetOpSplitP.__new__(_py_TetOpSplitP )
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_TetOpSplitP from_ref(const TetOpSplitP &ref):
        _py_TetOpSplitP.from_ptr(<TetOpSplitP*>&ref)
