####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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

try:
    import cPickle as Pickle
except ImportError:
    import Pickle

try:
    import steps
    HAS_STEPS = True
except ImportError:
    HAS_STEPS = False

ELEM_VERTEX = 0
ELEM_TRI = 1
ELEM_TET = 2
ELEM_UNDEFINED = 99

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

class ShadowComp:
    """
    Shadow class for STEPS compartment.
    
    Parameters:
        * name        Name of the compartment
        * mesh        ShadowMesh object
        * indices     Indices of the tetrahedrons in the compartment
        * vsyss       List of ids of the associated volume systems
    """
    def __init__(self, name, mesh, indices, vsyss):
        """
        Constructor.
        """
        self.name = name
        self.indices = indices
        if not isinstance(vsyss, list):
            raise Exception("vsyss should be a list of volume system ids.")
        self.vsyss = vsyss
        
        mesh.addComp(self)
    
    def exportTo(self, file):
        """
        Export data to a given file.
            
        Parameters:
            * file        File for the exporting
            
        Return:
            None
            """
        Pickle.dump(self.name, file)
        Pickle.dump(self.indices, file)
        Pickle.dump(self.vsyss, file)
    
    @staticmethod
    def importFrom(mesh, file):
        """
        Import data from a given file and create a ShadowComp object according to the data.
            
        Parameters:
            * mesh        Associated ShadowMesh file
            * file        File for the importing
            
        Return:
            ShadowComp object
        """
        name = Pickle.load(file)
        indices = Pickle.load(file)
        vsyss = Pickle.load(file)
        
        return ShadowComp(name, mesh, indices, vsyss)

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

class ShadowPatch:
    """
    Shadow class for STEPS patch.
    
    Parameters:
        * name        Name of the patch
        * mesh        ShadowMesh object
        * indices     Indices of the triangles in the patch
        * ssyss       List of ids of the associated surface systems
        * icomp       Inner ShadowComp object
        * ocomp       Outer ShadowComp object
    """
    def __init__(self, name, mesh, indices, ssyss, icomp, ocomp = None):
        """
        Constructor.
        
        """
        self.name = name
        self.indices = indices
        if not isinstance(ssyss, list):
            raise Exception("ssyss should be a list of surface system ids.")
        self.ssyss = ssyss
        self.icomp = icomp
        self.ocomp = ocomp
        
        mesh.addPatch(self)
    
    def exportTo(self, file):
        """
        Export data to a given file.
            
        Parameters:
            * file        File for the exporting
            
        Return:
            None
        """
        Pickle.dump(self.name, file)
        Pickle.dump(self.indices, file)
        Pickle.dump(self.ssyss, file)
        Pickle.dump(self.icomp.name, file)
        if self.ocomp != None:
            Pickle.dump(self.ocomp.name, file)
        else:
            Pickle.dump(None, file)
    
    @staticmethod
    def importFrom(mesh, file):
        """
        Import data from a given file and create a ShadowPatch object according to the data.
            
        Parameters:
            * mesh        Associated ShadowMesh file
            * file        File for the importing
            
        Return:
            ShadowPatch object
        """
        name = Pickle.load(file)
        indices = Pickle.load(file)
        ssyss = Pickle.load(file)
        icomp = mesh.comps[Pickle.load(file)]
        ocomp = Pickle.load(file)
        if ocomp != None:
            ocomp = mesh.comps[ocomp]
        return ShadowPatch(name, mesh, indices, ssyss, icomp, ocomp)


################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

class ShadowPatch:
    """
        Shadow class for STEPS patch.
        
        Parameters:
        * name        Name of the patch
        * mesh        ShadowMesh object
        * indices     Indices of the triangles in the patch
        * ssyss       List of ids of the associated surface systems
        * icomp       Inner ShadowComp object
        * ocomp       Outer ShadowComp object
        """
    def __init__(self, name, mesh, indices, ssyss, icomp, ocomp = None):
        """
            Constructor.
            
            """
        self.name = name
        self.indices = indices
        if not isinstance(ssyss, list):
            raise Exception("ssyss should be a list of surface system ids.")
        self.ssyss = ssyss
        self.icomp = icomp
        self.ocomp = ocomp
        
        mesh.addPatch(self)
    
    def exportTo(self, file):
        """
            Export data to a given file.
            
            Parameters:
            * file        File for the exporting
            
            Return:
            None
            """
        Pickle.dump(self.name, file)
        Pickle.dump(self.indices, file)
        Pickle.dump(self.ssyss, file)
        Pickle.dump(self.icomp.name, file)
        if self.ocomp != None:
            Pickle.dump(self.ocomp.name, file)
        else:
            Pickle.dump(None, file)
    
    @staticmethod
    def importFrom(mesh, file):
        """
            Import data from a given file and create a ShadowPatch object according to the data.
            
            Parameters:
            * mesh        Associated ShadowMesh file
            * file        File for the importing
            
            Return:
            ShadowPatch object
            """
        name = Pickle.load(file)
        indices = Pickle.load(file)
        ssyss = Pickle.load(file)
        icomp = mesh.comps[Pickle.load(file)]
        ocomp = Pickle.load(file)
        if ocomp != None:
            ocomp = mesh.comps[ocomp]
        return ShadowPatch(name, mesh, indices, ssyss, icomp, ocomp)


################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################
class ShadowMesh:
    """
        Shadow class for STEPS mesh.
        
        Parameters:
            None
    """
    
    def __init__(self):
        self.included_tets = set([])
        self.included_tris = set([])
        self.comps = {}
        self.patches = {}
        self.rois = {}
    
    def addComp(self, comp):
        """
        Add a compartment to the mesh.
        Note: Internal method, should use Comp(name, mesh, indices, vsys) instead.
        """
        
        if comp.name in self.comps:
            raise Exception("Compartment with name " + comp.name + " already exists.")
        
        intersect = self.included_tets & set(comp.indices)
        if len(intersect) != 0:
            raise Exception("The following tetrahedrons have been assigned to other compartment: " + str(intersect))
        
        self.comps[comp.name] = comp
        self.included_tets = self.included_tets | set(comp.indices)
    
    def addPatch(self, patch):
        """
        Add a patch to the mesh.
        Note: Internal method, should use Patch(name, mesh, indices, ssys, icomp, ocomp) instead.
        """
        
        if patch.name in self.patches:
            raise Exception("Patch with name " + patch.name + " already exists.")
        
        intersect = self.included_tris & set(patch.indices)
        if len(intersect) != 0:
            raise Exception("The following triangles have been assigned to other patch: " + str(intersect))
        
        if patch.icomp.name not in self.comps:
            raise Exception("Compartment " + patch.icomp.name + " does not exist.")
        
        if patch.ocomp != None:
            if patch.ocomp.name != None and (patch.ocomp.name not in self.comps):
                raise Exception("Compartment " + patch.ocomp.name + " does not exist.")
        
        self.patches[patch.name] = patch
        self.included_tris = self.included_tris | set(patch.indices)
    
    def addROI(self, roi_name, type, indices):
        """
        Add Region of Interest Data.
            
        Parameters:
            * roi_name        Name of the region of interest
            * type            Type of the elements
            * indices         Indices of the elements
            
        Return:
            None
        """
        if roi_name in self.rois:
            raise Exception("Region of Interest " + roi_name + " already exists.")
        self.rois[roi_name] = {"Type":type, "Indices":set(indices)}
    
    def removeROI(self, roi_name):
        """
        Remove Region of Interest Data.
            
        Parameters:
            * roi_name        Name of the region of interest
            
        Return:
            None
        """
        if roi_name not in self.rois:
            raise Exception("Region of Interest " + roi_name + " does not exist.")
        self.rois.pop(roi_name)
    
    def exportTo(self, file_name):
        """
        Export data to a given file.
            
        Parameters:
            * file        File for the exporting
            
        Return:
            None
        """
        file = open(file_name, "wb")
        Pickle.dump(len(self.comps), file)
        Pickle.dump(len(self.patches), file)
        for c in self.comps.values():
            c.exportTo(file)
        for p in self.patches.values():
            p.exportTo(file)
        Pickle.dump(self.rois, file)
        file.close()
    
    @staticmethod
    def importFrom(file_name):
        """
        Import data from a given file and create a ShadowMesh object according to the data.
            
        Parameters:
            * file        File for the importing
            
        Return:
            ShadowMesh object
        """
        file = open(file_name, "rb")
        mesh = ShadowMesh()
        ncomps = Pickle.load(file)
        npatches = Pickle.load(file)
        for c in range(ncomps):
            comp = ShadowComp.importFrom(mesh, file)
        for p in range(npatches):
            patch = ShadowPatch.importFrom(mesh, file)
        mesh.rois = Pickle.load(file)
        file.close()
        return mesh
    
    def writeToTetmesh(self, tetmesh, node_proxy, tet_proxy, tri_proxy):
        """
        Write data to STEPS Tetmesh object.
            
        Parameters:
            * tetmesh     STEPS Tetmesh object
            * node_proxy  ElementProxy object for nodes
            * tet_proxy   ElementProxy object for tetrahedrons
            * tri_proxy   ElementProxy object for triangles
            
        Return:
            ShadowMesh object
            
            """
        if not HAS_STEPS:
            raise Exception("This function is not available without STEPS integration with CUBIT.")
        else:
            temp_comps = {}
            for c in self.comps.values():
                print "Write compartment ", c.name, " to Tetmesh."
                steps_indices = [tetproxy.getSTEPSID(i) for i in c.indices]
                comp = sgeom.TmComp(c.name, tetmesh, steps_indices)
                temp_comps[c.name] = comp
                for v in c.vsyss:
                    comp.addVolsys(v)
            if len(self.patches) != 0:
                if tri_proxy.getSize() == 0:
                    print "ImportAbaqus does not support triangle index mapping, use ImportAbaqus2 instead."
                else:
                    for p in self.patches.values():
                        print "Write patch ", p.name, " to Tetmesh."
                        steps_indices = [triproxy.getSTEPSID(i) for i in c.indices]
                        icomp = temp_comps[p.icomp]
                        ocomp = None
                        if p.ocomp != None:
                            ocomp = temp_comps[p.ocomp]
                        patch = sgeom.TmPatch(p.name, tetmesh, steps_indices, icomp, ocomp)
                        for s in c.ssyss:
                            patch.addSurfsys(s)
            for roi in self.rois.items():
                roi_name = roi[0]
                roi_type = roi[1]["Type"]
                roi_import_indices = roi[1]["Indices"]
                print "Write ROI data ", roi[0], " to Tetmesh."
                if roi_type == ELEM_VERTEX:
                    steps_indices = [nodeproxy.getSTEPSID(i) for i in roi_import_indices]
                elif roi_type == ELEM_TET:
                    steps_indices = [tetproxy.getSTEPSID(i) for i in roi_import_indices]
                elif roi_type == ELEM_TRI:
                    if tri_proxy.getSize() == 0:
                        print "ImportAbaqus does not support triangle index mapping, use ImportAbaqus2 instead."
                    else:
                        steps_indices = [triproxy.getSTEPSID(i) for i in roi_import_indices]
                else:
                    steps_indices = roi_import_indices
                tetmesh.addROI(roi_name, roi_type, steps_indices)


################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################
