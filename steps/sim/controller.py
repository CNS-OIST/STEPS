# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class FuncCore(object):
    
    """This class implements the core functionality of a solver module.
    """

    
    def _rsf(self, funcname):
        """Resolve a function, specified by funcname, in the solver module."""
        return self._solver_module.__dict__.get(funcname, None)
    
    
    def __init__(self, solver_module, model, geom, rng):
        # Copy a refence to the solver module.
        self._solver_module = solver_module
        
        # Resolve core module functions of the solver interface.
        self._siNewState = self._rsf('siNewState')
        self._siDelState = self._rsf('siDelState')
        self._siBeginStateDef = self._rsf('siBeginStateDef')
        self._siEndStateDef = self._rsf('siEndStateDef')
        self._siSetRNG = self._rsf('siSetRNG')
        self._siBeginVarDef = self._rsf('siBeginVarDef')
        self._siEndVarDef = self._rsf('siEndVarDef')
        self._siNewSpec = self._rsf('siNewSpec')
        self._siBeginReacDef = self._rsf('siBeginReacDef')
        self._siEndReacDef = self._rsf('siEndReacDef')
        self._siNewReac = self._rsf('siNewReac')
        self._siAddReacLHS = self._rsf('siAddReacLHS')
        self._siAddReacRHS = self._rsf('siAddReacRHS')
        self._siBeginCompDef = self._rsf('siBeginCompDef')
        self._siEndCompDef = self._rsf('siEndCompDef')
        self._siNewComp = self._rsf('siNewComp')
        self._siAddCompSpec = self._rsf('siAddCompSpec')
        self._siAddCompReac = self._rsf('siAddCompReac')
        self._siReset = self._rsf('siReset')
        self._siRun = self._rsf('siRun')
        self._siGetTime = self._rsf('siGetTime')
        self._siGetCompVol = self._rsf('siGetCompVol')
        self._siSetCompVol = self._rsf('siSetCompVol')
        self._siGetCompCount = self._rsf('siGetCompCount')
        self._siSetCompCount = self._rsf('siSetCompCount')
        self._siGetCompMass = self._rsf('siGetCompMass')
        self._siSetCompMass = self._rsf('siSetCompMass')
        self._siGetCompConc = self._rsf('siGetCompConc')
        self._siSetCompConc = self._rsf('siSetCompConc')
        self._siGetCompClamped = self._rsf('siGetCompClamped')
        self._siSetCompClamped = self._rsf('siSetCompClamped')
        self._siGetCompReacKf = self._rsf('siGetCompReacKf')
        self._siSetCompReacKf = self._rsf('siSetCompReacKf')
        self._siGetCompReacActive = self._rsf('siGetCompReacActive')
        self._siSetCompReacActive = self._rsf('siSetCompReacActive')

        # Check whether the essential interface functions were successfully
        # resolved.
        assert self._siNewState != None
        assert self._siDelState != None
        # ... (TODO)

        # Now, attempt to create a state.
        self._state = self._siNewState()
        self._siSetRNG(self._state, rng)
        
        # Setup phase. Pass the model and geometric structure to the
        # the simulator, allow it to define a state.
        # TODO: re-organize this w.r.t. child class
        self._siBeginStateDef(self._state)
        
        # Setup species.
        self._lut_specs = { }
        self._setupVars(model)
        
        # Initialize: declare all reactions. 
        self._lut_reacs = { }
        self._setupReacs(model)
        
        # Initialize: declare all compartments.
        self._lut_comps = { }
        self._setupComps(model, geom)
        
        # Finish the state definition, create actual state(?)
        # Get ready for simulation.
        # TODO: again, reorganize (think about mesh and stuff)
        self._siEndStateDef(self._state)
    
    
    def __del__(self):
        self._siDelState(self._state)
        self._state = None

    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def _spec(self, spec):
        if type(spec) is str:
            spec = self._lut_specs[spec]
        return spec

    def _reac(self, reac):
        if type(reac) is str:
            reac = self._lut_reacs[reac]
        return reac

    def _comp(self, comp):
        if type(comp) is str:
            comp = self._lut_comps[comp]
        return comp

    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def _setupVars(self, model):
        self._siBeginVarDef(self._state)
        for s in model.getAllSpecies():
            sid = s.id
            sidx = self._siNewSpec(self._state, sid)
            self._lut_specs[sid] = sidx
        self._siEndVarDef(self._state)
	
    
    def _setupReacs(self, model):
        self._siBeginReacDef(self._state)
        for vsys in model.getAllVolsys():
            for reac in vsys.getAllReactions():
                rid = reac.id
                ridx = self._siNewReac(self._state, rid)
                self._lut_reacs[rid] = ridx
                for lhs in reac.lhs:
                    self._siAddReacLHS(self._state, \
                        ridx, self._spec(lhs.id))
                for rhs in reac.rhs:
                    self._siAddReacRHS(self._state, \
                        ridx, self._spec(rhs.id))
        self._siEndReacDef(self._state)
    
    
    def _setupComps(self, model, geom):
        self._siBeginCompDef(self._state)
        for c in geom.getAllComps():
            cid = c.id
            cidx = self._siNewComp(self._state, cid)
            self._lut_comps[cid] = cidx
            for vsys in c.volsys:
                # Find the Volsys object with the given name.
                vsys = model.getVolsys(vsys)
                for spec in vsys.getAllSpecies():
                    self._siAddCompSpec(self._state, \
                        cidx, self._spec(spec.id))
                for reac in vsys.getAllReactions():
                    self._siAddCompReac(self._state, \
                        cidx, self._reac(reac.id))
        self._siEndCompDef(self._state)
        

	#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #


    def reset(self):
        self._siReset(self._state)
		

    def run(self, endtime):
        self._siRun(self._state, endtime)


	#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #


    def getTime(self):
        return self._siGetTime(self._state)

    time = property(getTime)

		
	#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
	
    
    def getCompVol(self, comp):
        return self._siGetCompVol(self._state, self._comp(comp))
    
    def setCompVol(self, comp, vol):
        self._siSetCompVol(self._state, self._comp(comp), vol)
    
    
    def getCompCount(self, comp, spec):
        return self._siGetCompCount(self._state, \
            self._comp(comp), self._spec(spec))
    
    def setCompCount(self, comp, spec, num):
        self._siSetCompCount(self._state, \
            self._comp(comp), self._spec(spec), num)


    def getCompMass(self, comp, spec):
        return self._siGetCompMass(self._state, \
            self._comp(comp), self._spec(spec))
    
    def setCompMass(self, comp, spec, mass):
        self._siSetCompMass(self._state, \
            self._comp(comp), self._spec(spec), mass)
    

    def getCompConc(self, comp, spec):
        return self._siGetCompConc(self._state, \
            self._comp(comp), self._spec(spec))
    
    def setCompConc(self, comp, spec, conc):
        self._siSetCompConc(self._state, \
            self._comp(comp), self._spec(spec), conc)
    

    def getCompClamped(self, comp, spec):
        return self._siGetCompClamped(self._state, \
            self._comp(comp), self._spec(spec))
    
    def setCompClamped(self, comp, spec, buf):
        self._siGetCompClamped(self._state, \
            self._comp(comp), self._spec(spec), buf)


    def getCompReacKf(self, comp, reac):
        return self._siGetCompReacKf(self._state, \
            self._comp(comp), self._reac(reac))

    def setCompReacKf(self, comp, reac, kf):
        self._siSetCompReacKf(self._state, \
            self._comp(comp), self._reac(reac), kf)

    
    def getCompReacActive(self, comp, reac):
        return self._siGetCompReacActive(self._state, \
            self._comp(comp), self._reac(reac))
    
    def setCompReacActive(self, comp, reac, act):
        self._siSetCompReacActive(self._state, \
            self._comp(comp), self._reac(reac), act)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class FuncSSA(FuncCore):

    def __init__(self):
        # Resolve SSA-specific functionality in the solver module interface.
        self._siStep = self._rsf('siStep')
        self._siGetNSteps = self._rsf('siGetNSteps')
        self._siGetA0 = self._rsf('siGetA0')
        self._siGetCompReacFC = self._rsf('siGetCompReacFC')
        self._siGetCompReacFH = self._rsf('siGetCompReacFH')
        self._siGetCompReacFA = self._rsf('siGetCompReacFA')
        self._siGetCompReacExtent = self._rsf('siGetCompReacExtent')
        self._siResetCompReacExtent = self._rsf('siResetCompReacExtent')

    def step(self):
        return self._siStep(self._state)

    def getNSteps(self):
        return self._siGetNSteps(self._state)

    nsteps = property(getNSteps)
    
    def getA0(self):
        return self._siGetA0(self._state)
        
    a0 = property(getA0)
    
    def getCompReacFC(self, comp, reac):
        return self._siGetCompReacFC(self._state, \
            self._comp(comp), self._reac(reac))
    
    def getCompReacFH(self, comp, reac):
        return self._siGetCompReacFH(self._state, \
            self._comp(comp), self._reac(reac))
    
    def getCompReacFA(self, comp, reac):
        return slf._siGetCompReacFA(self._state, \
            self._comp(comp), self._reac(reac))
    
    def getCompReacExtent(self, comp, reac):
        return self._siGetCompReacExtent(self._state, \
            self._comp(comp), self._reac(reac))
    
    def resetCompReacExtent(self, comp, reac):
        self._siResetCompReacExtent(self._state, \
            self._comp(comp), self._reac(reac))
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    
# END
