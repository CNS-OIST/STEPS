####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

import numpy as np
import pyqtgraph as pg

import steps.API_1.visual as svisual
from steps.API_1.geom import INDEX_DTYPE

from . import utils as nutils
from . import model as nmodel
from . import geom as ngeom
from . import sim as nsim
from . import saving as nsaving

__all__ = [
    'SimControl',
    'PlotDisplay',
    'TimePlot',
    'SpatialPlot',
    'NewRow',
    'SimDisplay',
    'ElementDisplay',
    'PartitionDisplay',
]

###################################################################################################
# Exposed API


class SimControl(nutils.UsableObject):
    """Main visual class to wrap a simulation

    :param title: Title of the simulation control window
    :type title: str
    :param start_time: Initial start time of the simulation
    :type start_time: float
    :param end_time: Initial end time of the simulation
    :type end_time: float
    :param upd_interval: Initial update interval of the simulation
    :type upd_interval: float

    The SimControl object should be used as a context manager for declaring displays and displays
    should be used as context managers for declaring specific plots or 3D elements::

        sim = Simulation(...)

        rs = ResultSelector(sim)              # We need a result selector path root for describing
                                              # things to add to the plots.

        ...                                   # Setting up initial conditions

        sc = SimControl()                     # Creating the SimControl object

        with sc:                              # Then used as a context-manager

            with PlotDisplay('Species plot'): # Declaring a PlotDisplay in the sc SimControl and
                                              # using it as a context-manager to declare
                                              # sub-components.

                TimePlot(rs.comp1.S1.Count)   # First row will show a plot of the number of S1 in
                                              # comp1 as a function of time.

                NewRow()                      # Next plots will be added to a different row

                SpatialPlot(                  # Second row will show the distribution of S1 along
                    rs.TETS().S1.Count,       # the z axis.
                    axis=[0, 0, 1],
                    nbins=100
                )

            with SimDisplay('3D plots'):      # Declaring a SimDisplay in the sc SimControl and
                                              # using it as a context-manager to declare which
                                              # elements should appear in the SimDisplay.

                ElementDisplay(               # All compartment and patches should be displayed
                    rs.ALL(Comp, Patch)
                )
                ElementDisplay(               # The S1 species in all compartments should be plotted
                    rs.ALL(Comp).S1,          # in yellow
                    color=(0, 0.5, 0.5, 1)
                )
                ElementDisplay(               # The S2 species in all patches should be plotted in
                    rs.ALL(Patch).S2,         # purple.
                    color=(0.5, 0, 0.5, 1)
                )
                ElementDisplay(               # Complexes CC in comp1 should be plotted with a
                    rs.comp1.CC,              # color that depends on their state.
                    color=
                        lambda x: (
                            x.Count(C1),
                            x.Count(C2),
                            x.Count(C3),
                            1,
                        ),
                )

        sc.run()                              # Launch the simulation control windows
    """

    def __init__(self, title="Sim Control", start_time=0.0, end_time=10, upd_interval=0.1, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._app = pg.mkQApp()

        self.title = title
        self.start_time = start_time
        self.end_time = end_time
        self.upd_interval = upd_interval

        self._stepsSimControl = None

    def __exit__(self, exc_type, exc_val, exc_tb):
        super().__exit__(exc_type, exc_val, exc_tb)
        plotdisps = [pd._stepsPlotDisplay for pd in self._getChildrenOfType(PlotDisplay)]
        simdisps = [sd._stepsSimDisplay for sd in self._getChildrenOfType(SimDisplay)]
        allSolvers = [sim.solver for sim in self._getAllSims()]
        if len(allSolvers) > 0:
            self._stepsSimControl = svisual.SimControl(
                allSolvers, simdisps, plotdisps, self.title, self.start_time, self.end_time, self.upd_interval
            )

    def _getAllSims(self):
        sims = set()
        for c in self._getChildrenOfType(PlotDisplay, SimDisplay):
            sims |= c._getAllSims()
        return sims

    def run(self):
        """Launch the simulation displays"""
        self._window = self._app.exec_()

    def ALL(self, *cls):
        """:meta private:"""
        raise NotImplementedError()

    def __getattr__(self, name):
        """:meta private:"""
        raise AttributeError()


class PlotDisplay(nutils.UsingObjects(SimControl), nutils.UsableObject):
    """A window containing one or several plots

    :param title: Title of the window
    :type title: str
    :param size: Size of the window, in pixels
    :type size: Tuple[int, int]

    The PlotDisplay object should be used as a context manager for declaring specific plots. See
    example in :py:class:`SimControl`. Several plots can be declared in the same display. By
    default, the plots are added on a single row, from left to right. A new row can be added by
    calling :py:class:`NewRow`. Subsequent plots will be added to the new row.
    """

    def __init__(self, title=None, size=(800, 600), *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._stepsPlotDisplay = svisual.PlotDisplay(title, size)

    def _getAllSims(self):
        sims = set()
        for _, c in self.children.items():
            sims |= c._getAllSims()
        return sims

    def ALL(self, *cls):
        """:meta private:"""
        raise NotImplementedError()

    def __getattr__(self, name):
        """:meta private:"""
        raise AttributeError()


class TimePlot(nutils.UsingObjects(PlotDisplay)):
    """Plot time-dependent values

    :param rspath: Result selector path of the values to be plotted
    :type rspath: :py:class:`ResultSelector`
    :param title: Title of the plot
    :type title: str
    :param pen: Pen used to draw the plot lines
    :type pen: :py:class:`pyqtgraph.QPen`
    :param data_size: Maximum number of time points that should be displayed per line
    :type data_size: int
    :param x_range: Range of the x axis, automatic if None
    :type x_range: Union[None, Tuple[float, float]]
    :param y_range: Range of the y axis, automatic if None
    :type y_range: Union[None, Tuple[float, float]]
    :param show_x_grid: Display x axis grid
    :type show_x_grid: bool
    :param show_y_grid: Display y axis grid
    :type show_y_grid: bool
    :param x_label: Label for the x axis (time)
    :type x_label: str
    :param y_label: Label for the y axis, automatic if None (based on ``rspath`` labels)
    :type y_label: Union[None, str]
    :param label_style: Label style parameters
    :type label_style: dict
    :param \*\*kwargs: Other keywords supported by :py:class:`pyqtgraph.PlotItem`

    See example in :py:class:`SimControl`. The result selector path can contain several values.
    If no pen is specified, the different values will be automatically attributed different colors,
    otherwise, all lines will share the same pen.
    """

    def __init__(
        self,
        rspath,
        title=None,
        pen=None,
        data_size=1000,
        x_range=None,
        y_range=None,
        show_x_grid=True,
        show_y_grid=True,
        label_pen=None,
        x_label=('Time', 's'),
        y_label=None,
        label_style={'color': '#ffffff', 'font-size': '16px'},
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        if not isinstance(rspath, nsaving.ResultSelector):
            raise TypeError(f'Expected a ResultSelector, got {rspath} instead.')
        if title is None:
            title = rspath._strDescr()

        (disp,) = self._getUsedObjects()
        sdisp = disp._stepsPlotDisplay

        self._rspath = rspath

        if title in sdisp.updater:
            raise NameError(f'A plot with name {title} already exists.')
        plot = sdisp.addPlot(title=title)
        sdisp.updater[title] = _GenericTimePlotUpdater(
            plot, rspath, data_size, x_range, y_range, pen, **kwargs
        )

        lblpen = pg.mkPen(color=(255, 255, 255), width=2) if label_pen is None else label_pen
        plot.getAxis('left').setPen(lblpen)
        plot.getAxis('bottom').setPen(lblpen)
        plot.showGrid(x=show_x_grid, y=show_y_grid)
        if x_label is not None:
            plot.setLabel('bottom', *x_label, **label_style)
        if y_label is not None:
            plot.setLabel('left', *y_label, **label_style)

        self._stepsPlot = plot

    def _getAllSims(self):
        return set([self._rspath.sim])

    def ALL(self, *cls):
        """:meta private:"""
        raise NotImplementedError()

    def __getattr__(self, name):
        """:meta private:"""
        raise AttributeError()


class SpatialPlot(nutils.UsingObjects(PlotDisplay)):
    """Plot space-dependent values

    :param rspath: Result selector path of the values to be plotted
    :type rspath: :py:class:`ResultSelector`
    :param title: Title of the plot
    :type title: str
    :param mode: ``'distr'`` or ``'mean'``, see explanation below
    :type mode: str
    :param axis: Axis on which the data should be projected
    :type axis: Tuple[float, float, float]
    :param nbins: Number of bins
    :type nbins: int
    :param x_label: Label for the x axis, first element is the label title, second is the unit
    :type x_label: Tuple[str, str]
    :param y_range: Range of the y axis, automatic if None
    :type y_range: Union[None, Tuple[float, float]]
    :param show_x_grid: Display x axis grid
    :type show_x_grid: bool
    :param show_y_grid: Display y axis grid
    :type show_y_grid: bool
    :param label_pen: Pen used to draw the labels
    :type label_pen: :py:class:`pyqtgraph.QPen`
    :param label_style: Label style parameters
    :type label_style: dict
    :param \*\*kwargs: Other keywords supported by :py:class:`pyqtgraph.PlotItem`

    See example in :py:class:`SimControl`. The result selector path should encompass values in
    different places in space (i.e using ``TETS(...)`` or ``TRIS()``). The position of each of the
    covered geometrical elements is then projected on the given axis and data is binned with the
    given number of bins. If ``'distr'`` mode is used, the plot will be a histogram displaying the
    accumulated (summed) values in each bin. If ``'mean'`` mode is used, the plot will display the
    average value for each bin.
    """

    def __init__(
        self,
        rspath,
        title=None,
        mode='distr',
        axis=np.array([1, 0, 0]),
        nbins=20,
        x_label=('Position', 'm'),
        y_range=None,
        show_x_grid=True,
        show_y_grid=True,
        label_pen=None,
        label_style={'color': '#ffffff', 'font-size': '16px'},
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        locTypes = rspath.metaData['loc_type']
        locIds = rspath.metaData['loc_id']
        mesh = rspath.sim.geom
        if all(tpe == ngeom.TetReference._locStr for tpe in locTypes):
            elemLst = ngeom.TetList(locIds, mesh=mesh)
        elif all(tpe == ngeom.TriReference._locStr for tpe in locTypes):
            elemLst = ngeom.TriList(locIds, mesh=mesh)
        else:
            raise TypeError(
                f'Cannot create a SpatialPlot from a result selector that does '
                f'not explicitely correspond to tetrahedrons or triangles.'
            )

        if title is None:
            title = rspath._strDescr()

        (disp,) = self._getUsedObjects()
        sdisp = disp._stepsPlotDisplay

        self._rspath = rspath

        if title in sdisp.updater:
            raise NameError(f'A plot with name {title} already exists.')
        plot = sdisp.addPlot(title=title)
        sdisp.updater[title] = _GenericSpatialPlotUpdater(
            plot, rspath, nbins, axis, elemLst, y_range, mode, **kwargs
        )

        pen = pg.mkPen(color=(255, 255, 255), width=2) if label_pen is None else label_pen
        plot.getAxis('left').setPen(pen)
        plot.getAxis('bottom').setPen(pen)
        plot.showGrid(x=show_x_grid, y=show_y_grid)
        if x_label is not None:
            plot.setLabel('bottom', *x_label, **label_style)

        self._stepsPlot = plot

    def _getAllSims(self):
        return set([self._rspath.sim])

    def ALL(self, *cls):
        """:meta private:"""
        raise NotImplementedError()

    def __getattr__(self, name):
        """:meta private:"""
        raise AttributeError()


class NewRow(nutils.UsingObjects(PlotDisplay)):
    """Add a new row to the display"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        (disp,) = self._getUsedObjects()
        disp._stepsPlotDisplay.nextRow()

    def _getAllSims(self):
        return set()

    def ALL(self, *cls):
        """:meta private:"""
        raise NotImplementedError()

    def __getattr__(self, name):
        """:meta private:"""
        raise AttributeError()


class SimDisplay(nutils.UsingObjects(SimControl), nutils.UsableObject):
    """A window containing 3D elements

    :param title: Title of the window
    :type title: str
    :param size: Size of the window, in pixels
    :type size: Tuple[int, int]

    The SimDisplay object should be used as a context manager for declaring specific elements. See
    example in :py:class:`SimControl`.
    """

    def __init__(self, title=None, size=(800, 600), *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._stepsSimDisplay = svisual.SimDisplay(title, w=size[0], h=size[1])

    def add(self, elemDisp):
        """Add a previously declared element display

        :param elemDisp: The element display
        :type elemDisp: :py:class:`ElementDisplay`

        This is useful when wanting to reuse element displays that were declared in other
        SimDisplays.
        """
        self._addChildren(elemDisp)
        for sve in elemDisp._stepsVisElems:
            self._stepsSimDisplay.addItem(sve)

    def merge(self, *simDisp):
        """Merge with a previously declared simulation display

        :param \*simDisp: The simulation displays to merge to self
        :type \*simDisp: :py:class:`SimDisplay`

        This method adds all :py:class:`ElementDisplay` that were associated to the given
        simulation display(s) to the current simulation display. This is useful when separate
        views of subcomponents are necessary but the user still wants a full view using all
        components.
        """
        for sd in simDisp:
            for ed in sd._getChildrenOfType(ElementDisplay):
                self.add(ed)

    def _getAllSims(self):
        sims = set()
        for _, c in self.children.items():
            sims |= c._getAllSims()
        return sims

    def ALL(self, *cls):
        """:meta private:"""
        raise NotImplementedError()

    def __getattr__(self, name):
        """:meta private:"""
        raise AttributeError()


class ElementDisplay(nutils.UsingObjects(SimDisplay)):
    """3D element to be displayed in a SimDisplay

    :param path: A result selector path to a geometrical location or to an object (see
        :py:class:`steps.API_2.sim.SimPath`)
    :type path: :py:class:`steps.API_2.saving.ResultSelector`
    :param color: A color for the element or a function that takes an element as a parameter and
        returns a color (in case several elements are provided). An automatic color is selected if
        this is None.
    :type color: Union[Callable[[StepsElement], Tuple[float, ...]], Tuple[float, ...], None]
    :param spec_size: If species are provided as elements, defines their size in the 3D plot.
    :type spec_size: float
    :param max_nspec: Maximum number of species that can be visualized in the plot.
    :type max_nspec: int
    :param max_density: Maximum density of species that can be visualized in the plot.
    :type max_density: float
    :param auto_adjust: Flag for auto adjustment of visualized species counts
    :type auto_adjust: bool

    See example in :py:class:`SimControl`. If the result selector path that is given points to
    geometrical elements (like compartments, patches, tetrahedrons, etc.), the corresponding 3D
    surface is added to the display. If the result selector path points to objects instead (like
    species, complexes, channels, etc.), points corresponding to their location are added to the
    display.
    """

    def __init__(
        self,
        path,
        color=None,
        spec_size=0.2,
        max_nspec=10000,
        max_density=1000e18,
        auto_adjust=True,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        (disp,) = self._getUsedObjects()
        sdisp = disp._stepsSimDisplay

        if isinstance(path, nsaving._ResultPath):
            path = path.simpath
        if not isinstance(path, nsim.SimPath):
            raise TypeError(f'Expected a SimPath, got {path} instead.')

        self._sim = path._sim

        self._stepsVisElems = []
        for *_loc, elem in path._walk(expand=False, combine=False):
            tris = ElementDisplay._getSurfTrisFromElem(elem)
            # First try to use the element as a volume
            if tris is not None:
                self._stepsVisElems.append(_GenericSurfElem(elem, tris, sdisp, color, *args, **kwargs))
            elif isinstance(
                elem, (nmodel.Species, nmodel.Complex, nmodel.ComplexSelector, nmodel.ComplexState)
            ):
                # Transform location in list of tets or tris
                *_, loc = _loc
                tets = ElementDisplay._getTetsFromElem(loc)
                tris = ElementDisplay._getTrisFromElem(loc) if tets is None else None
                for subElem in elem._simPathWalkExpand():
                    if tets is not None:
                        self._stepsVisElems.append(
                            _GenericTetScatterElem(
                                tets,
                                (*_loc, subElem),
                                sdisp,
                                color,
                                spec_size,
                                max_nspec,
                                max_density,
                                auto_adjust,
                                *args,
                                **kwargs,
                            )
                        )
                    elif tris is not None:
                        self._stepsVisElems.append(
                            _GenericTriScatterElem(
                                tris,
                                (*_loc, subElem),
                                sdisp,
                                color,
                                spec_size,
                                max_nspec,
                                max_density,
                                auto_adjust,
                                *args,
                                **kwargs,
                            )
                        )
                    else:
                        raise TypeError(f'Could not extract spatial information from {loc}.')

    def _getAllSims(self):
        return set([self._sim])

    @classmethod
    def _getSurfTrisFromElem(cls, elem):
        if isinstance(elem, ngeom.TetReference):
            return elem.faces
        elif isinstance(elem, ngeom.TriReference):
            return ngeom.TriList([elem])
        elif isinstance(elem, ngeom.TetList):
            return elem.surface
        elif isinstance(elem, ngeom.TriList):
            return elem
        elif hasattr(elem, 'tets'):
            return elem.tets.surface
        elif hasattr(elem, 'tris'):
            return elem.tris
        return None

    @classmethod
    def _getTetsFromElem(cls, elem):
        if isinstance(elem, ngeom.TetReference):
            return ngeom.TetList([elem])
        elif isinstance(elem, ngeom.TetList):
            return elem
        elif hasattr(elem, 'tets'):
            return elem.tets
        return None

    @classmethod
    def _getTrisFromElem(cls, elem):
        if isinstance(elem, ngeom.TriReference):
            return ngeom.TriList([elem])
        elif isinstance(elem, ngeom.TriList):
            return elem
        elif hasattr(elem, 'tris'):
            return elem.tris
        return None

    def ALL(self, *cls):
        """:meta private:"""
        raise NotImplementedError()

    def __getattr__(self, name):
        """:meta private:"""
        raise AttributeError()


class PartitionDisplay(nutils.UsingObjects(SimControl)):
    """A window displaying a mesh partition

    :param partition: The partition to be displayed
    :type partition: :py:class:`steps.API_2.geom.MeshPartition`
    :param elem: ``'tet'`` for displaying tetrahedron partition and ``'tri'`` for displaying
        triangle partition
    :type elem: str
    :param title: Title of the window
    :type title: str
    :param size: Size of the window, in pixels
    :type size: Tuple[int, int]
    """

    def __init__(self, partition, elem='tet', title=None, size=(800, 600), *args, **kwargs):
        super().__init__(*args, **kwargs)
        if elem == 'tet':
            self._stepsPartitionDisplay = svisual.TetPartitionDisplay(
                partition._mesh.stepsMesh, partition._tet_hosts, title, w=size[0], h=size[1], *args, **kwargs
            )
        elif elem == 'tri':
            self._stepsPartitionDisplay = svisual.TriPartitionDisplay(
                partition._mesh.stepsMesh, partition._tri_hosts, title, w=size[0], h=size[1], *args, **kwargs
            )
        else:
            raise ValueError(f'Parameter elem needs to be equal to "tet" or "tri", got "{elem}" instead.')

    def _getAllSims(self):
        return set()

    def ALL(self, *cls):
        """:meta private:"""
        raise NotImplementedError()

    def __getattr__(self, name):
        """:meta private:"""
        raise AttributeError()


###################################################################################################
# Internal generic plots


class _GenericTimePlotUpdater:
    """
    Generic plot updater, replaces all the plot updaters from visual/Plotting.py
    """

    COLOR_NB = 5
    COLOR_IND = 0

    def __init__(self, plot, rspath, data_size, x_range, y_range, pen, **kwargs):
        self.plot = plot
        self.plotkwargs = kwargs
        self.rspath = rspath
        self.data_size = data_size

        if rspath._getEvalLen() > 1:
            self.plot.addLegend()

        self.time = [rspath.sim.Time]
        self.data = [[v] for v in rspath._evaluate()]
        self.curves = []
        for dat, label in zip(self.data, self.rspath.labels):
            if pen is None:
                npen = pg.intColor(
                    _GenericTimePlotUpdater.COLOR_IND,
                    hues=_GenericTimePlotUpdater.COLOR_NB,
                )
                _GenericTimePlotUpdater.COLOR_IND += 1
            else:
                npen = pen
            self.curves.append(plot.plot(self.time, dat, name=label, pen=npen, **kwargs))

        if x_range is not None:
            self.plot.setXRange(x_range[0], x_range[1])
        if y_range is not None:
            self.plot.setYRange(y_range[0], y_range[1])

    def update(self):
        self.time.append(self.rspath.sim.Time)
        for i, v in enumerate(self.rspath._evaluate()):
            self.data[i].append(v)
        if len(self.time) > self.data_size:
            self.time.pop(0)
            for dat in self.data:
                dat.pop(0)
        for curve, dat in zip(self.curves, self.data):
            curve.setData(self.time, dat)

    def reset(self):
        self.time = [self.rspath.sim.Time]
        self.data = [[v] for v in self.rspath._evaluate()]
        for curve, dat in zip(self.curves, self.data):
            curve.setData(self.time, dat)


class _GenericSpatialPlotUpdater:
    """
    Generic plot updater, replaces all the dist updaters from visual/Plotting.py
    """

    def __init__(self, plot, rspath, nbins, axis, elems, y_range, mode='distr', **kwargs):
        self.plot = plot
        self.rspath = rspath
        self.mode = mode

        positions = [elem.center @ axis for elem in elems]
        bins = np.linspace(min(positions), max(positions), nbins + 1)
        self.dig = np.digitize(positions, bins)
        # Move the max values to the last bin
        self.dig[self.dig == max(self.dig)] -= 1
        # Shift indices to discard bin number 0 (lower than min)
        self.dig = np.array([v - 1 for v in self.dig])
        y_data = self._getYData()

        self.x_data = [bins[0]]
        for b in bins[1:-1]:
            self.x_data += [b] * 2
        self.x_data.append(bins[-1])

        if self.mode == 'distr':
            self.curve = plot.plot(
                self.x_data, y_data, fillLevel=0, fillBrush=kwargs.get('fillBrush', (0, 0, 128)), **kwargs
            )
        elif self.mode == 'mean':
            self.curve = plot.plot(self.x_data, y_data, **kwargs)
        else:
            raise NotImplementedError(f'Mode {self.mode} is not implemented.')

        if y_range is not None:
            self.plot.setYRange(y_range[0], y_range[1])

    def update(self):
        self.curve.setData(self.x_data, self._getYData())

    def reset(self):
        self.update()

    def _doubleValues(self, data):
        res = []
        for v in data:
            res += [v] * 2
        return res

    def _getYData(self):
        if self.mode == 'distr':
            return self._doubleValues(np.bincount(self.dig, weights=self.rspath._evaluate()))
        elif self.mode == 'mean':
            bincounts = np.bincount(self.dig)
            data = np.bincount(self.dig, weights=self.rspath._evaluate())
            data[bincounts > 0] = data[bincounts > 0] / bincounts[bincounts > 0]
            data[bincounts == 0] = 0
            return self._doubleValues(data)
        else:
            raise NotImplementedError(f'Mode {self.mode} is not implemented.')


class _GenericSurfElem(pg.opengl.GLMeshItem):
    DEFAULT_ALPHA = 0.1
    COLOR_NB = 5
    COLOR_IND = 0

    def __init__(self, elem, tris, sdisp, color=None, *args, **kwargs):
        self.id = elem.name
        self.display = sdisp

        if color is None:
            ncolor = pg.intColor(
                _GenericSurfElem.COLOR_IND,
                hues=_GenericSurfElem.COLOR_NB,
                alpha=_GenericSurfElem.DEFAULT_ALPHA * 255,
            )
            _GenericSurfElem.COLOR_IND += 1
        else:
            ncolor = color

        self.bound_min = tris.mesh.bbox.min * self.display.scale
        self.bound_max = tris.mesh.bbox.max * self.display.scale

        steps_mesh = tris.mesh._getStepsObjects()[0]
        surface_tris = np.array([tri.idx for tri in tris], dtype=INDEX_DTYPE)
        v_set_size = steps_mesh.getTriVerticesSetSizeNP(surface_tris)
        tris_data = np.zeros(surface_tris.size * 3, dtype=INDEX_DTYPE)
        v_set = np.zeros(v_set_size, dtype=INDEX_DTYPE)
        verts_data = np.zeros(v_set_size * 3)
        steps_mesh.getTriVerticesMappingSetNP(surface_tris, tris_data, v_set)
        steps_mesh.getBatchVerticesNP(v_set, verts_data)
        verts_data *= self.display.scale
        tris_data.shape = -1, 3
        verts_data.shape = -1, 3

        mesh_data = pg.opengl.MeshData(vertexes=verts_data, faces=tris_data)
        pg.opengl.GLMeshItem.__init__(
            self,
            meshdata=mesh_data,
            smooth=False,
            computeNormals=True,
            shader='balloon',
            glOptions='additive',
        )
        self.setColor(ncolor)
        self.display.addItem(self)

    def updateItem(self):
        return


class _GenericScatterElem(pg.opengl.GLScatterPlotItem):
    DEFAULT_ALPHA = 1
    COLOR_NB = 10
    COLOR_IND = 0

    def __init__(
        self, elems, path, sdisp, color, spec_size, max_nspec, max_density, auto_adjust, *args, **kwargs
    ):
        *_loc, elem = path

        self.id = elem.name
        self.display = sdisp
        self.smesh = elems.mesh._getStepsObjects()[0]

        if len(path) != 3:
            raise Exception(f'Invalid path for plotting elements in tetrahedrons: {path}')

        self.bound_min = elems.mesh.bbox.min * self.display.scale
        self.bound_max = elems.mesh.bbox.max * self.display.scale

        if color is None:
            self.col = pg.colorTuple(
                pg.intColor(
                    _GenericScatterElem.COLOR_IND,
                    hues=_GenericScatterElem.COLOR_NB,
                    alpha=_GenericScatterElem.DEFAULT_ALPHA * 255,
                )
            )
            _GenericScatterElem.COLOR_IND += 1
        elif hasattr(color, '__call__'):
            self.col = color(elem)
        else:
            self.col = color

        self.spec_size = spec_size
        self.max_nspec = max_nspec
        self.max_density = max_density
        self.auto_adjust = auto_adjust

        self.elems = np.array([tet.idx for tet in elems], dtype=INDEX_DTYPE)
        self.rspath = self._getResSelect(path)

        data = self._getData()

        pg.opengl.GLScatterPlotItem.__init__(
            self, pos=data, color=self.col, size=self.spec_size, pxMode=False
        )
        self.display.addItem(self)

    def _getData(self):
        points = np.array(self.rspath._evaluate(), dtype=np.uint32)
        total = np.sum(points)
        if total > self.max_nspec:
            total = self.__reduce(points)
        data = np.zeros(3 * int(total))
        self._genElemVisualPointsNP(self.elems, points, data)
        data *= self.display.scale
        data.shape = -1, 3
        return data

    def __reduce(self, point_counts):
        """
        Reduce the number of points being generated.
        """
        self._reduceBatchPointCountsNP(self.elems, point_counts, self.max_density)
        total = np.sum(point_counts)
        temp_density = self.max_density
        while total > self.max_nspec and self.auto_adjust:
            temp_density *= self.max_nspec / total
            self._reduceBatchPointCountsNP(self.elems, point_counts, temp_density)
            total = np.sum(point_counts)
        return total

    def updateItem(self):
        """
        Update the component.
        """
        self.setData(pos=self._getData(), color=self.col)

    def _getResSelect(self, path):
        pass

    def _genElemVisualPointsNP(self, *args):
        pass

    def _reduceBatchPointCountsNP(self, *args):
        pass


class _GenericTetScatterElem(_GenericScatterElem):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _getResSelect(self, path):
        return nsaving.ResultSelector(path[0]).TETS(self.elems).LIST(path[2]).Count

    def _genElemVisualPointsNP(self, *args):
        self.smesh.genTetVisualPointsNP(*args)

    def _reduceBatchPointCountsNP(self, *args):
        self.smesh.reduceBatchTetPointCountsNP(*args)


class _GenericTriScatterElem(_GenericScatterElem):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _getResSelect(self, path):
        return nsaving.ResultSelector(path[0]).TRIS(self.elems).LIST(path[2]).Count

    def _genElemVisualPointsNP(self, *args):
        self.smesh.genTriVisualPointsNP(*args)

    def _reduceBatchPointCountsNP(self, *args):
        self.smesh.reduceBatchTriPointCountsNP(*args)
