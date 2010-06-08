
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2009 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#  Last Changed Rev:  $Rev: 322 $
#  Last Changed Date: $Date: 2010-04-07 10:43:03 +0900 (Wed, 07 Apr 2010) $
#  Last Changed By:   $Author: wchen $

""" 
.. Note::

    This module is preliminary, means some of the functions are still under development.
    Code modification / debugging is wellcomed.
    Please email steps.dev@gmail.com if you would like to share you changes with others.

Visual Toolkit

"""

try:
    import wx
    import sys
    import  wx.lib.colourselect as  csel
    haveWXPython = True
except ImportError:
    haveWXPython = False


try:
    from wx import glcanvas
    haveGLCanvas = True
except ImportError:
    haveGLCanvas = False

try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *
    haveOpenGL = True
except ImportError:
    haveOpenGL = False

try:
    import steps
    from steps import geom
    haveSTEPS = True
except ImportError:
    haveSTEPS = False

import time
import os   
import cPickle 

wildcard = "STEPS CheckPoint File (*.checkpoint)|*.checkpoint|"     \
           "All files (*.*)|*.*"
#----------------------------------------------------------------------
#
#                          Mesh Display
#
#----------------------------------------------------------------------

class MyCanvasBase(glcanvas.GLCanvas):
    def __init__(self, parent):
        glcanvas.GLCanvas.__init__(self, parent, -1)
        self.init = False
        # initial mouse position
        self.lastx = self.x = 0
        self.lasty = self.y = 0
        self.rotate = [20.0, 45.0, 0.0]
        self.centerx = 0.0
        self.centery = 0.0
        self.zooming = 1.0
        self.size = None
        self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.Bind(wx.EVT_LEFT_DOWN, self.OnMouseDown)
        self.Bind(wx.EVT_LEFT_UP, self.OnMouseUp)
        self.Bind(wx.EVT_MOTION, self.OnMouseMotion)
        self.Bind(wx.EVT_MOUSEWHEEL, self.OnMouseWheel)
        
    def OnEraseBackground(self, event):
        pass # Do nothing, to avoid flashing on MSW.


    def OnSize(self, event):
        size = self.size = self.GetClientSize()
        if self.GetContext():
            self.SetCurrent()
            glViewport(0, 0, size.width, size.height)
        event.Skip()


    def OnPaint(self, event):
        dc = wx.PaintDC(self)
        self.SetCurrent()
        if not self.init:
            self.InitGL()
            self.init = True
        self.OnDraw()
        
    def OnMouseDown(self, evt):
        self.freezexy = False
        self.CaptureMouse()
        self.x, self.y = self.lastx, self.lasty = evt.GetPosition()


    def OnMouseUp(self, evt):
        self.ReleaseMouse()
    
    def OnMouseMotion(self, evt):
        if evt.Dragging() and evt.LeftIsDown():
            self.lastx, self.lasty = self.x, self.y
            self.x, self.y = evt.GetPosition()
            if self.x > self.lastx:
                self.centerx += 0.5 * self.zooming
            elif self.x < self.lastx:
                self.centerx -= 0.5 * self.zooming
            if self.y > self.lasty:
                self.centery -= 0.5 * self.zooming
            elif self.y < self.lasty:
                self.centery += 0.5 * self.zooming
            self.Refresh(False)
                
    def OnMouseWheel(self, evt):
        if evt.GetWheelRotation() > 0:
            self.zooming = self.zooming * 0.9
        else:
            self.zooming = self.zooming * 1.1
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(self.boundmin * self.zooming, self.boundmax * self.zooming, self.boundmin * self.zooming, self.boundmax * self.zooming, self.boundmin, self.boundmax)
        glMatrixMode (GL_MODELVIEW)
        self.OnDraw()

    def ResetDisp(self):
        self.rotate = [0.0, 0.0, 0.0]
        self.centerx = 0.0
        self.centery = 0.0
        self.zooming = 1.0
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(self.boundmin * self.zooming, self.boundmax * self.zooming, self.boundmin * self.zooming, self.boundmax * self.zooming, self.boundmin, self.boundmax)
        glMatrixMode (GL_MODELVIEW)
        self.OnDraw()
    
class MeshCanvas(MyCanvasBase):
    def __init__(self, parent, mesh, solver, scale = 1e-6):
        self.parent = parent
        self.solver = solver
        self.mesh = mesh
        self.scale = scale
        self.bgcolor = [255, 255, 255]
        self.meshcolor = [0, 0, 0]
        self.meshalpha = 255
        self.specs_data = {}
        self.pointsize = 4.0
        self.__loadMesh__()
        self.updateSpecies(parent.spec_color_mapping)
        MyCanvasBase.__init__(self, parent)

    def __loadMesh__(self):
        self.ntets = self.mesh.ntets
        self.nverts = self.mesh.nverts
        self.ntris = self.mesh.ntris
        self.verts = []
        self.tris = []
        self.tribounds = self.mesh.getTriBoundary()
        boundmin = self.mesh.getBoundMin()
        boundmax = self.mesh.getBoundMax()
        bmin = min(boundmin)
        bmax = max(boundmax)
        self.boundmin = bmin * 1.2 / self.scale
        self.boundmax = bmax * 1.2 / self.scale

        for n in range(self.nverts):
            unscaled = self.mesh.getVertex(n)
            self.verts.append([unscaled[0]/self.scale,unscaled[1]/self.scale,unscaled[2]/self.scale])
        for n in range(self.ntris):
            self.tris.append(self.mesh.getTri(n))

    def InitGL(self):
        # set viewing projection
        glMatrixMode(GL_PROJECTION)
        glOrtho(self.boundmin, self.boundmax, self.boundmin, self.boundmax, self.boundmin, self.boundmax)

        # position viewer
        glMatrixMode(GL_MODELVIEW)
        #glTranslatef(0.0, 0.0, -2.0)

        glEnable(GL_DEPTH_TEST)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glEnable( GL_BLEND )
        #glEnable(GL_LIGHTING)
        #glEnable(GL_LIGHT0)


    def OnDraw(self):
        # clear color and depth buffers
        glClearColor(self.bgcolor[0]/255.0, self.bgcolor[1]/255.0, self.bgcolor[2]/255.0, 1.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()      
        
        glTranslatef(self.centerx, self.centery, 0)
        
        glRotatef(self.rotate[0], 1.0, 0.0, 0.0)
        glRotatef(self.rotate[1], 0.0, 1.0, 0.0) 
        glRotatef(self.rotate[2], 0.0, 0.0, 1.0) 
        
        # draw species
        glEnableClientState(GL_VERTEX_ARRAY)
        glPointSize(self.pointsize)
        for key in self.specs_data.keys():
            glColor3ubv(self.parent.spec_color_mapping[key][1])
            glVertexPointer(3, GL_DOUBLE, 0, self.specs_data[key])
            nverts = len(self.specs_data[key]) / 3
            glDrawElements(GL_POINTS, nverts, GL_UNSIGNED_BYTE, range(nverts)) 
               
        glDisableClientState(GL_VERTEX_ARRAY)
        
        glPointSize(1)        
        # draw mesh
        colorvec = [self.meshcolor[0],self.meshcolor[1],self.meshcolor[2],self.meshalpha]
        glColor4ubv(colorvec)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
        glBegin(GL_TRIANGLES)
        for n in range(len(self.tribounds)):
            tri = self.tris[self.tribounds[n]]
            for v in range(3):
                vert = tri[v]
                glVertex3f( self.verts[vert][0], self.verts[vert][1], self.verts[vert][2])
        glEnd()

        glFlush()
        self.SwapBuffers()
    
    def updateSpecies(self, spec_color_mapping):
        self.specs_data.clear()
        for key in spec_color_mapping.keys():
            if spec_color_mapping[key][0]:
                self.specs_data[key] = []
                
        for t in range(self.ntets):
            for key in self.specs_data.keys():
                if(self.solver.getTetCount(t, key) > 0):
                    center = self.mesh.getTetBarycenter(t)
                    scaled_center = [center[0]/self.scale, center[1]/self.scale, center[2]/self.scale]
                    self.specs_data[key].extend(scaled_center)
        

#----------------------------------------------------------------------
#
#                           GUI Simulation Frontend
#
#----------------------------------------------------------------------

class SimCtrlPanel(wx.Panel):
    def __init__(self, parent, model, mesh, solver):
        wx.Panel.__init__(self, parent)
        self.parent = parent
        self.model = model
        self.mesh = mesh
        self.solver = solver
        self.spec_list = []
        
        # layout
        box = wx.BoxSizer(wx.VERTICAL)
        
        box.AddSpacer(10)
        box.Add(wx.StaticText(self, -1, "Simulation Time (Sec)"), 0, wx.ALIGN_CENTER, border = 10)
        self.simtime_txt = wx.StaticText(self, -1, "%e" % (self.solver.getTime()))
        box.AddSpacer(10)
        box.Add(self.simtime_txt, 0, wx.ALIGN_CENTER, border = 10)
                      
        static1 = wx.StaticText(self, -1, "Stop at (Sec)")
        self.stoptime_text = wx.TextCtrl(self,-1, "1e-3")
        box.AddSpacer(10)
        box.Add(static1, 0, wx.ALIGN_CENTER, border = 10)
        box.AddSpacer(10)
        box.Add(self.stoptime_text, 0, wx.ALIGN_CENTER, border = 10)
        
        
        static2 = wx.StaticText(self, -1, "Advance Interval (Sec)")
        self.advance_text = wx.TextCtrl(self,-1, "1e-5")
        
        box.AddSpacer(10)
        box.Add(static2, 0, wx.ALIGN_CENTER, border = 10)
        box.AddSpacer(10)
        box.Add(self.advance_text, 0, wx.ALIGN_CENTER, border = 10)
        

        self.runbtn = wx.Button(self, -1, "Run")
        box.AddSpacer(10)
        box.Add(self.runbtn, 0, wx.ALIGN_CENTER, border = 10)
        
        #self.replaybtn = wx.Button(self, -1, "Replay")
        #box.AddSpacer(10)
        #box.Add(self.replaybtn, 0, wx.ALIGN_CENTER, border = 10)

        box.AddSpacer(10)
        btn_box = color_box = wx.BoxSizer(wx.HORIZONTAL)
        self.savebtn = wx.Button(self, -1, "Save")
        btn_box.AddSpacer(5)
        btn_box.Add(self.savebtn, 0, wx.ALIGN_CENTER, border = 10)
        
        self.loadbtn = wx.Button(self, -1, "Load")
        btn_box.AddSpacer(5)
        btn_box.Add(self.loadbtn, 0, wx.ALIGN_CENTER, border = 10)
        
        box.Add(btn_box, 0, wx.ALIGN_CENTER, border = 5)
        
        box.AddSpacer(20)
        box.Add(wx.StaticText(self, -1, "Species"), 0, wx.ALIGN_CENTER, border = 10)
        box.AddSpacer(10)
        self.spec_listbox = wx.ListBox(self, -1, style = wx.LB_SINGLE)
        box.Add(self.spec_listbox, 1, wx.ALIGN_CENTER, border = 10)
        
        
        color_box = wx.BoxSizer(wx.HORIZONTAL)
        color_box.AddSpacer(10)
        self.specdisp_check = wx.CheckBox(self, -1, "Display")
        color_box.Add(self.specdisp_check, 0, wx.ALIGN_CENTER, border = 10)
        self.speccolor_select = csel.ColourSelect(self,-1, colour = wx.RED)
        color_box.AddSpacer(10)
        color_box.Add(self.speccolor_select, 0, wx.ALIGN_CENTER, border = 10)
        self.specsize_spin = wx.SpinCtrl(self, -1, "Size")
        self.specsize_spin.SetRange(1,100)
        self.specsize_spin.SetValue(4)
        
        
        box.AddSpacer(10)
        box.Add(color_box, 0, wx.ALIGN_CENTER, border = 10)
        
        size_box = wx.BoxSizer(wx.HORIZONTAL)
        size_box.AddSpacer(10)
        size_box.Add(wx.StaticText(self, -1, "Spec Size"), 0, wx.ALIGN_CENTER, border = 10)
        size_box.AddSpacer(10)
        size_box.Add(self.specsize_spin, 0, wx.ALIGN_CENTER, border = 10)
        
        box.AddSpacer(5)
        box.Add(size_box, 0, wx.ALIGN_CENTER, border = 10)

        mesh_box = wx.BoxSizer(wx.HORIZONTAL)
        mesh_box.Add(wx.StaticText(self, -1, "Mesh Color"), 0, wx.ALIGN_CENTER, border = 10)
        self.meshcolor_select = csel.ColourSelect(self,-1, colour = wx.BLACK)
        mesh_box.AddSpacer(20)
        mesh_box.Add(self.meshcolor_select, 0, wx.ALIGN_CENTER, border = 10)
        
        box.AddSpacer(20)
        box.Add(mesh_box, 0, wx.ALIGN_CENTER, border = 10)
        
        box.AddSpacer(5)
        box.Add(wx.StaticText(self, -1, "Transparency"), 0, wx.ALIGN_CENTER, border = 10)
        
        self.alpha_slider = wx.Slider(self, -1, 255, 0, 255, style = wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
        box.Add(self.alpha_slider, 0, wx.ALIGN_CENTER, border = 10)
        
        background_box = wx.BoxSizer(wx.HORIZONTAL)
        background_box.Add(wx.StaticText(self, -1, "BG Color"), 0, wx.ALIGN_CENTER, border = 10)
        self.backgroundcolor_select = csel.ColourSelect(self,-1, colour = wx.WHITE)
        background_box.AddSpacer(10)
        background_box.Add(self.backgroundcolor_select, 0, wx.ALIGN_CENTER, border = 10)
        
        box.AddSpacer(10)
        box.Add(background_box, 0, wx.ALIGN_CENTER, border = 10)
        
        box.AddSpacer(10)
        box.Add(wx.StaticText(self, -1, "Rotation"), 0, wx.ALIGN_CENTER, border = 10)
        self.rotatebox = wx.RadioBox(self, -1, "", wx.DefaultPosition, wx.DefaultSize,["X","Y", "Z"], 1, wx.RA_SPECIFY_ROWS)
        box.Add(self.rotatebox, 0, wx.ALIGN_CENTER, border = 10)
        self.rotate_slider = wx.Slider(self, -1, 20, -180, 180, style = wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
        box.Add(self.rotate_slider, 0, wx.ALIGN_CENTER, border = 10)
        self.resetbtn = wx.Button(self, -1, "Reset")
        box.Add(self.resetbtn, 0, wx.ALIGN_CENTER, border = 10)
        
        box.AddSpacer(20)             
        self.SetSizerAndFit(box)
        
        # event binding
        self.meshcolor_select.Bind(csel.EVT_COLOURSELECT, self.OnMeshColorChange, id=self.meshcolor_select.GetId())
        self.alpha_slider.Bind(wx.EVT_SCROLL, self.OnMeshAlphaChange, id = self.alpha_slider.GetId())
        self.backgroundcolor_select.Bind(csel.EVT_COLOURSELECT, self.OnBGColorChange, id=self.backgroundcolor_select.GetId())
        self.spec_listbox.Bind(wx.EVT_LISTBOX, self.OnSpecSelectChange, id = self.spec_listbox.GetId())
        self.specdisp_check.Bind(wx.EVT_CHECKBOX, self.OnSpecDispChange, id = self.specdisp_check.GetId())
        self.speccolor_select.Bind(csel.EVT_COLOURSELECT, self.OnSpecColorChange, id=self.speccolor_select.GetId())
        self.specsize_spin.Bind(wx.EVT_SPINCTRL, self.OnSpecSizeChange, id = self.specsize_spin.GetId())
        self.runbtn.Bind(wx.EVT_BUTTON, self.OnRun, id = self.runbtn.GetId())
        self.rotatebox.Bind(wx.EVT_RADIOBOX, self.OnCordChange, id = self.rotatebox.GetId())
        self.rotate_slider.Bind(wx.EVT_SCROLL, self.OnRotateDisp, id = self.rotate_slider.GetId()) 
        self.resetbtn.Bind(wx.EVT_BUTTON, self.OnResetView, id = self.resetbtn.GetId())
        self.savebtn.Bind(wx.EVT_BUTTON, self.OnSave, id = self.savebtn.GetId())
        self.loadbtn.Bind(wx.EVT_BUTTON, self.OnLoad, id = self.loadbtn.GetId())       
        #self.replaybtn.Bind(wx.EVT_BUTTON, self.OnReplay, id = self.replaybtn.GetId()) 
        
    def refreshSpecList(self, list):
        self.spec_list = list
        self.spec_listbox.Set(list)
    
    def OnMeshColorChange(self, event):
        self.parent.disp_panel.meshcolor = list(event.GetValue())
        self.parent.disp_panel.OnDraw()
        
    def OnMeshAlphaChange(self, event):
        self.parent.disp_panel.meshalpha = self.alpha_slider.GetValue()
        self.parent.disp_panel.OnDraw()
    
    def OnBGColorChange(self, event):
        self.parent.disp_panel.bgcolor = list(event.GetValue())
        self.parent.disp_panel.OnDraw()
        
    def OnSpecSelectChange(self, event):
        specstatus = self.parent.spec_color_mapping[event.GetString()]
        self.specdisp_check.SetValue(specstatus[0])
        self.speccolor_select.SetColour(specstatus[1])
    
    def OnSpecDispChange(self, event):
        specid = self.spec_listbox.GetSelections()
        if len(specid) == 1:
            specid = specid[0]
        else:
            return
        if specid >= 0 and specid < len(self.spec_list):
            self.parent.spec_color_mapping[self.spec_list[specid]][0] = event.IsChecked()
            self.parent.updateSpecs()
            self.parent.disp_panel.OnDraw()
        
    def OnSpecColorChange(self, event):
        specid = self.spec_listbox.GetSelections()
        if len(specid) == 1:
            specid = specid[0]
        else:
            return
        if specid >= 0 and specid < len(self.spec_list):
            self.parent.spec_color_mapping[self.spec_list[specid]][1] = list(event.GetValue())
            self.parent.updateSpecs()
            self.parent.disp_panel.OnDraw()
        
    def OnSpecSizeChange(self,event):
        self.parent.disp_panel.pointsize = self.specsize_spin.GetValue()
        self.parent.disp_panel.OnDraw()
    
    def OnCordChange(self, event):
        cord = event.GetInt()
        self.rotate_slider.SetValue(self.parent.disp_panel.rotate[cord])
    
    def OnRotateDisp(self, event):
        cord = self.rotatebox.GetSelection()
        self.parent.disp_panel.rotate[cord] = self.rotate_slider.GetValue()
        self.parent.disp_panel.OnDraw()
    
    def OnResetView(self, event):
        self.parent.disp_panel.ResetDisp()
        self.rotate_slider.SetValue(0)
    
    def OnRun(self, event):
        while self.solver.getTime() < float(self.stoptime_text.GetValue()):
            self.solver.advance(float(self.advance_text.GetValue()))
            self.simtime_txt.SetLabel("%e" % (self.solver.getTime()))
            self.parent.updateSpecs()
            self.parent.disp_panel.OnDraw()
            self.Update()
            
    def OnSave(self,event):
        dlg = wx.FileDialog(
            self, message="Save file as ...", defaultDir=os.getcwd(), 
            defaultFile="Untitled.checkpoint", wildcard=wildcard, style=wx.SAVE
            )

        dlg.SetFilterIndex(2)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.solver.checkpoint(path)
        dlg.Destroy()
        
    def OnLoad(self, event):
        dlg = wx.FileDialog(
            self, message="Choose a file",
            defaultDir=os.getcwd(), 
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.CHANGE_DIR
            )

        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.solver.restore(path)
            self.simtime_txt.SetLabel("%e" % (self.solver.getTime()))
            self.parent.updateSpecs()
            self.parent.disp_panel.OnDraw()
        dlg.Destroy()
 
    #def OnReplay(self, event):
    #    dlg = wx.FileDialog(
    #        self, message="Choose a file",
    #        defaultDir=os.getcwd(), 
    #        defaultFile="",
    #        style=wx.OPEN | wx.CHANGE_DIR
    #        )
#
 #       if dlg.ShowModal() == wx.ID_OK:
  #          path = dlg.GetPath()
    #        self.solver.restore(path)
   #         self.simtime_txt.SetLabel("%e" % (self.solver.getTime()))
    #        self.parent.updateSpecs()
    #        self.parent.disp_panel.OnDraw()
     #   dlg.Destroy()
        
class SimPanel(wx.SplitterWindow):
    def __init__(self, parent, model, mesh, solver, scale = 1e-6, color_map = [[255,0,0], [0,0,255]]):
        wx.SplitterWindow.__init__(self, parent, -1)
        
        #get species list and initialise color mapping
        species = model.getAllSpecs()
        self.spec_color_mapping = {}
        
        for s in range(len(species)):
            if s < len(color_map):
                self.spec_color_mapping[species[s].id] = [True, color_map[s]]
            else:
                self.spec_color_mapping[species[s].id] = [True, [0,255,0]]
                
        self.SetAutoLayout(True)
        
        self.ctrl_panel = SimCtrlPanel(self, model, mesh, solver)
        self.ctrl_panel.refreshSpecList(self.spec_color_mapping.keys())
        self.disp_panel = MeshCanvas(self, mesh, solver, scale)
        
        self.SetMinimumPaneSize(100)
        self.SplitVertically(self.disp_panel, self.ctrl_panel, -100)
        
    def updateSpecs(self):
        self.disp_panel.updateSpecies(self.spec_color_mapping)
    

def GUISim(model, mesh, solver, scale = 1e-6, color_map = [[255,0,0], [0,0,255]]):
    """
    Create a graphical frontend for a mesh based simulation.
    
    .. Note::
    
        The graphical frontend should be used for mesh based simulation only.
    
    Arguements:
        * steps.model.Model model
        * steps.geom.Tetmesh mesh
        * steps.solver.Tetexact solver
        
    Example:
        see examples/diffusion
    """
    
    if not haveSTEPS:
        print("Cannot find STEPS! Please install it first.")
        return
    
    if not haveWXPython:
        print("Cannot find WxPython! Please install it first.")
        return
        
    if not haveOpenGL:
        print("Cannot find PyOpenGL! Please install it first.")
        return
        
        
    app = wx.App(False)
    frame = wx.Frame(None, wx.ID_ANY, 'STEPS Visual Frontend for Mesh Based Simulation http://steps.sourceforge.net', size = (1200,900), style=wx.DEFAULT_FRAME_STYLE)
    sim_panel = SimPanel(frame, model, mesh, solver, scale, color_map)
    frame.Centre()
    frame.Show()
    app.MainLoop()
    