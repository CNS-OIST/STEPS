####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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

import sys, time
from PyQt4 import QtCore, QtGui
import os
import sys

from threading import Thread


class WorkThread(QtCore.QThread):
    workDoneSignal = QtCore.pyqtSignal()
    def __init__(self, sim, endtime):
        QtCore.QThread.__init__(self)
    
        self.sim = sim
        self.endtime = endtime
    
    def __del__(self):
        self.wait()
    
    def run(self):
        self.sim.run(self.endtime)
        self.workDoneSignal.emit()
        return

class SimControl(QtGui.QWidget):
    """
    Simulation Control for the visualization toolkit
            
    Parameters:
        * sims                List of STEPS simulation solvers
        * sim_displays        List of SimDisplay objects
        * plot_displays       List of PlotDisplay objects
        * title               Title of the control interface
        * start_time          Initial start time of the simulation
        * end_time            Initial end time of the simulation
        * upd_interval        Initial update interval of the simulation
    """
    
    def __init__(self, sims, sim_displays, plot_displays, title = "Sim Control", start_time = 0.0, end_time = 10, upd_interval = 0.1):
        """
        Constructor.

        """
        QtGui.QWidget.__init__(self)
        self.setWindowTitle(title)
        self.resize(300, 50)
        self.setWindowTitle(title)
        self.sims = sims
        self.nsims = len(sims)
        self.sim_displays = sim_displays
        self.plot_displays = plot_displays
        visual_items = []
        for d in sim_displays:
            visual_items.extend(d.getItems())
        self.visual_items = set(visual_items)
        self.layout = QtGui.QGridLayout(self)
        self.workers = []
        
        self.runButton = QtGui.QPushButton("run")
        self.runButton.released.connect(self.__run)
        self.stopButton = QtGui.QPushButton("stop")
        self.stopButton.released.connect(self.__stop)
        self.stopButton.setDisabled(True)
        self.simForward = True
        self.checkpoint = False
        
        self.endtimeLabel = QtGui.QLabel("End Time(s): ")
        self.runtimeEdit = QtGui.QLineEdit(self)
        self.runtimeEdit.setText(str(end_time))
        self.updatetimeLabel = QtGui.QLabel("Update Interval(s): ")
        self.updateEdit = QtGui.QLineEdit(self)
        self.updateEdit.setText(str(upd_interval))
        self.simtimeLabel = QtGui.QLabel("Current Simulation Time: %fs" % (start_time))

        self.resetButton = QtGui.QPushButton("reset")
        self.resetButton.released.connect(self.__reset)
        
        self.__runSims(start_time, True, True)

        self.layout.addWidget(self.runButton, 0, 0)
        self.layout.addWidget(self.endtimeLabel, 0, 1)
        self.layout.addWidget(self.runtimeEdit, 0, 2)
        
        self.layout.addWidget(self.stopButton, 1, 0)
        self.layout.addWidget(self.updatetimeLabel, 1, 1)
        self.layout.addWidget(self.updateEdit, 1, 2)
        self.layout.addWidget(self.resetButton, 2, 0)
        self.layout.addWidget(self.simtimeLabel, 2, 1, 1, -1)
        self.show()

    def __runSims(self, endtime, once = False, checkpoint = False):
        self.finish_count = 0
        self.runButton.setDisabled(True)
        self.resetButton.setDisabled(True)
        self.stopButton.setEnabled(True)
        self.simForward = True
        for sim in self.sims:
            worker = WorkThread(sim, endtime)
            worker.workDoneSignal.connect(self.__workDone)
            self.workers.append(worker)
            worker.start()
        if once:
            self.simForward = False
        self.checkpoint = checkpoint

    def __run(self):
        end_time = float(self.runtimeEdit.text())
        update_interval = float(self.updateEdit.text())
        stop_time = self.current_time + update_interval
        
        if stop_time <= end_time and self.simForward:
            self.__runSims(stop_time)
        else:
            self.runButton.setEnabled(True)
            self.resetButton.setEnabled(True)
            self.stopButton.setDisabled(True)
            self.simForward = True

    
    def __workDone(self):
        self.finish_count += 1
        if self.nsims == self.finish_count:
            for worker in self.workers:
                worker.workDoneSignal.disconnect(self.__workDone)
            for item in self.visual_items:
                item.updateItem()
            for display in self.sim_displays:
                display.refresh()
            for display in self.plot_displays:
                display.refresh()
            self.current_time = self.sims[0].getTime()
            self.simtimeLabel.setText("Current Simulation Time: %.4es" % (self.current_time))
            self.workers = []
            if self.checkpoint:
                for i in range(len(self.sims)):
                    self.sims[i].checkpoint("_temp_sim%i.cp" % (i))
            self.__run()
    
    
    def __stop(self):
        self.runButton.setEnabled(True)
        self.resetButton.setEnabled(True)
        self.stopButton.setDisabled(True)
        self.simForward = False

    def __reset(self):
        for i in range(len(self.sims)):
            self.sims[i].restore("_temp_sim%i.cp" % (i))
        self.current_time = self.sims[0].getTime()
        self.simtimeLabel.setText("Current Simulation Time: %.4es" % (self.current_time))
        for item in self.visual_items:
            item.updateItem()
        for display in self.sim_displays:
            display.refresh()
        for display in self.plot_displays:
            display.reset()
            display.refresh()

    def addDisplay(self, d):
        self.sim_displays.append(d)

    def getsim_displays(self):
        return self.sim_displays

    def removeDisplay(self, d):
        self.sim_displays.remove(d)

    def hideDisplay(self, d):
        d.hide()

    def showDisplay(self, d):
        d.show()

    def closeEvent(self, event):
        for display in self.plot_displays:
            display.clear()
            display.hide()
            display.close()
        for display in self.sim_displays:
            display.close()
        
        event.accept()
  
