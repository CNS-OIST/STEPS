# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


"""Core classes and routines for describing model geometry.

A complete structural description of a model can be made with the classes 
in this module. Such descriptions can be used with well-mixed simulation
algorithms.
"""


import steps.error as serr
import steps.tools as stools


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Container(object):


    """
    """


    def __init__(self):
        """
        """
        self.__comps = { }
        self.__patches = { }


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _checkCompID(self, id):
        """Check if a given id is valid and not yet used by another 
        compartment.
        
        Parameters:
            id
                THe id that should be checked (must be a string).
        
        Returns:
            The id itself.
        
        Raises:
            steps.error.ArgumentError
                If the ID is not valid, or not unique.
        """
        stools.checkID(id)
        if id in self.__comps:
            raise serror.ArgumentError, '\'%s\' is already in use.' % id


    def _checkPatchID(self, id):
        """Check if a given id is valid and not yet used by another patch.
        
        Parameters:
            id
                THe id that should be checked (must be a string).
        
        Returns:
            The id itself.
        
        Raises:
            steps.error.ArgumentError
                If the ID is not valid, or not unique.
        """
        stools.checkID(id)
        if id in self.__patches:
            raise serror.ArgumentError, '\'%s\' is already in use.' % id


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleCompIDChange(self, oldid, newid):
        """Internal function.
        """
        if oldid == newid: return
        self._checkCompID(newid)
        c = self.__comps.pop(oldid)
        self.__comps[newid] = c


    def _handleCompAdd(self, comp):
        """Internal function.
        """
        assert comp._container == None, \
            '\'%s\' already assigned to a geometry container object.' % comp.id
        self._checkCompID(comp.id)
        comp._container = self
        self.__comps[comp.id] = comp


    def delComp(self, comp):
        """Delete a compartment from the geometry container.
        """
        comp = self.getComp(comp)
        self.__comps.pop(comp.id)
        comp._handleSelfDelete()


    def getComp(self, comp):
        """Resolve a compartment in the context of this geometric model.
        
        This method can be used to retrieve a compartment object by its id, 
        or to check whether a given compartment object belongs to this model.
        
        Arguments:
            comp
                Can be a string or a steps.geom.Comp object. If it is a
                string, the method will attempt to find the matching 
                compartment object. If it is a steps.geom.Comp object,
                the method will see if the object is part of this model.
            
        Returns:
            A steps.geom.Comp object.
        
        Raises:
            steps.error.ArgumentError
                If the id cannot be resolved, or if the compartment object
                does not belong to this model.
        """
        if isinstance(comp, str):
            try:
                comp = self.__comps[comp]
            except KeyError:
                raise serr.ArgumentError, \
                    'No compartment with id \'%s\'' % comp
        if comp._container != self:
            raise steps.ArgumentError, \
                'Compartment is no part of this geometry container.'
        return comp

        
    def getAllComps(self):
        """Return all compartments in this geometric model.
        """
        return self.__comps.values()


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
    
    
    def _handlePatchIDChange(self, oldid, newid):
        """Internal function.
        """
        if oldid == newid: return
        self._checkPatchID(newid)
        c = self.__patches.pop(oldid)
        self.__patches[newid] = c


    def _handlePatchAdd(self, patch):
        """Internal function.
        """
        assert patch._container == None, \
            '\'%s\' already assigned to a geometry container object.' %patch.id
        self._checkPatchID(patch.id)
        patch._container = self
        self.__patches[patch.id] = patch


    def delPatch(self, patch):
        """Delete a patch from the geometry container.
        """
        patch = self.getPatch(patch)
        self.__patches.pop(patch.id)
        patch._handleSelfDelete()


    def getPatch(self, patch):
        """Resolve a patch in the context of this geometric model.
        
        This method can be used to retrieve a patch object by its id, 
        or to check whether a given patch object belongs to this model.
        
        Arguments:
            patch
                Can be a string or a steps.geom.Patch object. If it is a
                string, the method will attempt to find the matching 
                patch object. If it is a steps.geom.Patch object,
                the method will see if the object is part of this model.
            
        Returns:
            A steps.geom.Patch object.
        
        Raises:
            steps.error.ArgumentError
                If the id cannot be resolved, or if the patch object does
                not belong to this model.
        """
        if isinstance(patch, str):
            try:
                patch = self.__patches[patch]
            except KeyError:
                raise serr.ArgumentError, \
                    'No patch with id \'%s\'' % patch
        if patch._container != self:
            raise steps.ArgumentError, \
                'Patch is no part of this geometry.'
        return patch

        
    def getAllPatches(self):
        """Return all patches in this geometric model.
        """
        return self.__patches.values()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Comp(object):


    """Base class for compartment objects.
    
    It provides basic functionality and data that is shared by all classes 
    derived from Comp: 
        * Getting and setting a valid compartment ID string, and handling
          the interaction with the container object.
        * Getting (and at least in this base class also setting) the total
          volume of the compartment.
        * The volume systems implemented in the compartment.
        * References to Patch objects.
    
    This base class can be used directly with well-mixed solvers.
    """


    def __init__(self, id, container, **params):
        """Initialize a compartment object.
        
        Parameters:
            id
                The ID of the compartment.
            container
                A reference to the container object.
                
            vol
                The volume of the compartment (might get ignored,
                depending on child class).
            volsys
                A sequence of volume systems.
        """
        self._container = None
        self.__id = stools.checkID(id)
        try:
            container._handleCompAdd(self)
        except:
            self._container = None
            raise
        assert self._container != None, \
            'Compartment not assigned to geometry container.'

        self._volsys = set()
        if 'volsys' in params: self.addVolsys(params['volsys'])
        self._vol = 0.0
        if 'vol' in params: self.vol = params['vol']
        
        self._ipatches = set()
        self._opatches = set()


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleSelfDelete(self):
        """Internal function.
        """
        self._container = None


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getID(self):
        """Return the compartment ID.
        """
        return self.__id

    def setID(self, id):
        """Change the compartment ID.
        """
        assert self._container != None, \
            'Compartment not assigned to geometry container.'
        if id == self.__id: return
        self._container._handleCompIDChange(self.__id, id)
        self.__id = id

    id = property(getID, setID)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getContainer(self):
        """Return the geometry container object.
        """
        return self._container

    container = property(getContainer)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def addVolsys(self, volsys):
        """Add a volume system to the compartment.
        
        Parameters:
            volsys
                A string or steps.model.Volsys object, or some sequence
                of such objects.
            
        Returns:
            ---
        
        Raises:
            ---
        """
        # Turn the input thingy into a set.
        volsys = set(volsys)
        # Now check its contents and add it.
        for vs in volsys:
            try:
                vs = vs.id
            except: pass
            self.volsys.add(vs)
    
    
    def delVolsys(self, volsys):
        """Remove a volume system from the compartment.
        
        Parameters:
            volsys
                A string or steps.model.Volsys object, or some sequence of
                such objects.
            
        Returns:
            ---
        
        Raises:
            ---
        """
        # Turn the input thingy into a set.
        volsys = set(volsys)
        # Now check its contents, and remove it.
        for vs in volsys:
            try:
                vs = vs.id
            except: pass
            self.volsys.discard(vs)


    def getVolsys(self):
        """Return a copy of the set of volume systems.
        
        By returning a (shallow) copy, we prevent meddling with this set.
        """
        return self._volsys.copy()
    
    volsys = property(getVolsys)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getVol(self):
        """
        """
        return self._vol

    def setVol(self, vol):
        """
        """
        # Python's float() type is actually a double..
        vol = float(vol)
        if vol <= 0.0:
            raise serr.ArgumentError, \
                'Compartment volume must be larger than zero.'
        self._vol = vol

    vol = property(getVol, setVol)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
    
    
    def _addIPatch(self, patch):
        """Internal function. 
        
        (Called from steps.geom.core.Patch.setOComp)
        """
        assert patch.container == self.container
        assert patch.ocomp == self
        self._ipatches.add(patch)
        
    
    def _delIPatch(self, patch):
        """Internal function.
        
        (Called from steps.geom.core.Patch.setOComp)
        """
        self._ipatches.discard(patch)
    
    
    def getIPatches(self):
        """Return a copy of the set of internal patches.
        
        By returning a (shallow) copy, we prevent meddling with this set.
        """
        return self._ipatches.copy()
    
    ipatches = property(getIPatches)
        
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _addOPatch(self, patch):
        """Internal function.
        
        (Called from steps.geom.core.Patch.setIComp)
        """
        assert patch.container == self.container
        assert patch.icomp == self
        self._opatches.add(patch)
    
    def _delOPatch(self, patch):
        """Internal function.
        
        (Called from steps.geom.core.Patch.setIComp)
        """
        self._opatches.discard(patch)
    
    
    def getOPatches(self):
        """Return a copy of the set of external patches.
        
        By returning a (shallow) copy, we prevent meddling with this set.
        """
        return self._opatches.copy()
    
    opatches = property(getOPatches)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Patch(object):


    """Base class for patch objects.
    
    A patch is a piece of 2D surface surrounding (part of) a 3D compartment.
    This base class provides basic functionality and descriptive data that 
    is shared by all types of patches ('type' meaning different types of 
    geometric descriptions):
    
        * Getting and setting a valid patch ID string, and handling
          the interaction with the container object.
          
        * Getting (and at least in this base class also setting) the total
          area of the patch.
          
        * The surface systems associated with the patches.
        
        * References to inside/outside compartments.
    
    This base class can be used directly with well-mixed solvers.
    """
    
    
    def __init__(self, id, container, icomp, ocomp = None):
        """
        """
        self._container = None
        self.__id = stools.checkID(id)
        try:
            container._handlePatchAdd(self)
        except:
            self._container = None
            raise
        assert self._container != None, \
            'Patch not assigned to geometry container.'

        # Needs to be re-written... what if it fails???
        self._icomp = None
        self.icomp = icomp
        self._ocomp = None
        self.ocomp = ocomp

        # Parse optional arguments.
        self._surfsys = set()
        if 'surfsys' in params: self.addSurfsys(params['surfsys'])
        self._area = 0.0
        if 'area' in params: self.area = params['area']


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleSelfDelete(self):
        """Internal function.
        """
        self._container = None


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getID(self):
        """Return the patch ID.
        """
        return self.__id

    def setID(self, id):
        """Change the patch ID.
        """
        assert self._container != None, \
            'Patch not assigned to geometry container.'
        if id == self.__id: return
        self._container._handlePatchIDChange(self.__id, id)
        self.__id = id

    id = property(getID, setID)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getContainer(self):
        """Return the container object.
        """
        return self._container

    container = property(getContainer)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def addSurfsys(self, surfsys):
        """Add surface system(s) to the patch.
        
        Parameters:
            surfsys
                A string or steps.model.Surfsys object, or some sequence
                of such objects.
            
        Returns:
            ---
        
        Raises:
            ---
        """
        # Turn the input thingy into a set.
        surfsys = set(surfsys)
        # Now check its contents and add it.
        for ss in surfsys:
            try:
                ss = ss.id
            except: pass
            self._surfsys.add(ss)
    
    
    def delSurfsys(self, surfsys):
        """Remove surface system(s) from the patch.
        
        Parameters:
            surfsys
                A string or steps.model.Surfsys object, or some sequence of
                such objects.
            
        Returns:
            ---
        
        Raises:
            ---
        """
        # Turn the input thingy into a set.
        surfsys = set(surfsys)
        # Now check its contents, and remove it.
        for ss in surfsys:
            try:
                ss = ss.id
            except: pass
            self._surfsys.discard(ss)
    
            
    def getSurfsys(self):
        """Return a copy of the set of surface systems.
        
        By returning a (shallow) copy, we prevent meddling with this set.
        """
        return self._surfsys.copy()
    
    surfsys = property(getSurfsys)
        
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getArea(self):
        """Return the surface area of the patch.
        """
        return self._area
    
    def setArea(self, area):
        """Set the surface area of the patch.
        
        May be overriden by a derived class, if this class compute the area 
        of the patch automatically (for instance by adding the areas of the
        constituent triangles).
        """
        area = float(area)
        if area <= 0.0:
            raise serr.ArgumentError, 'Patch area must be larger than zero.'
        self._area = area
    
    area = property(getArea, setArea)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
    
    
    def getIComp(self):
        """Return a reference to the inside compartment.
        """
        return self._icomp
    
    
    def _setIComp(self, icomp):
        """Set the inside compartment.
        
        Raises:
            steps.error.ArgumentError
        
        Notes:
            Should only be called during setup of the patch. Changing the
            internal compartment of a patch should not be changed at a later
            point in time.
        """
        # Do some tests on the specified compartment.
        if icomp.container != self.container:
            raise serr.ArgumentError, \
                'Compartment does not belong to same container as patch.'
        if self in icomp.ipatches:
            raise serr.ArgumentError, \
                'Patch is already on inside of compartment.'
        # Remove the patch if it was already on the outside of some other
        # compartment.
        if self._icomp != None:
            self._icomp._delOPatch(self)
        # Add it.
        self._icomp = icomp
        self._icomp._addOPatch(self)
    
    
    icomp = property(getIComp)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
    
    
    def getOComp(self):
        """Return a reference to the outside compartment.
        
        If the patch is a boundary patch, the method returns None.
        """
        return self._ocomp
    
    
    def _setOComp(self, ocomp):
        """Set the outside compartment.
        
        Raises:
            steps.error.ArgumentError
        
        Notes:
            Should only be called during setup of the patch. Changing the
            outer compartment of a patch should not be changed at a later
            point in time.
        """
        # Do some tests on the specified compartment.
        if ocomp.container != self.container:
            raise serr.ArgumentError, \
                'Compartment does not belong to same container as patch.'
        if self in ocomp.opatches:
            raise serr.ArgumentError, \
                'Patch is already on outside of compartment.'
        # Remove the patch if it was already on the inside of some other
        # compartment.
        if self._ocomp != None:
            self._ocomp._delIPatch(self)
        # Add it.
        self._ocomp = ocomp
        self._ocomp._addIPatch(self)
    
    
    ocomp = property(getOComp)
    

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
