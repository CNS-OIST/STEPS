####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

import argparse
import os.path
import pickle
import re
import subprocess
import tempfile
import typing

from . import utils
from stepsblender.blenderloader import HDF5BlenderLoader

SKIPPED_PARAMS = set([r'_', r'parent$', r'nameInParent$'])
DEFAULT_METAVARS = {
    int: 'i',
    float: 'x',
    str: 'string',
}


def addArgumentToParser(parser, fullName, tpe):
    kwargs = {}
    if typing.get_origin(tpe) is typing.Annotated:
        tpe, *annotations = typing.get_args(tpe)
        for annot in annotations:
            if isinstance(annot, dict):
                for key, val in annot.items():
                    if key in kwargs:
                        kwargs[key] += ' ' + val
                    else:
                        kwargs[key] = val
            elif 'help' in kwargs:
                kwargs['help'] += ' ' + annot
            else:
                kwargs['help'] = annot
    if tpe is bool:
        parser.add_argument(f'--{fullName}', action='store_true', default=argparse.SUPPRESS, **kwargs)
    else:
        kwargs.setdefault('metavar', DEFAULT_METAVARS.get(tpe, 'x'))
        parser.add_argument(f'--{fullName}', action='store', type=tpe, default=argparse.SUPPRESS, **kwargs)


def addVisualizationParameters(parser, params, prefix=tuple()):
    flatParams = []
    for name, param in params.items():
        if all(re.match(p, name) is None for p in SKIPPED_PARAMS):
            fullName = '.'.join(prefix + (name, ))
            if isinstance(param, dict):
                flatParams += addVisualizationParameters(parser, param, prefix + (name, ))
            else:
                defVal, tpe = param
                addArgumentToParser(parser, fullName, tpe)
                flatParams.append(fullName)
    return flatParams


def buildParameterDict(args, params):
    dct = {}
    for param in params:
        if param in args:
            splitParam = param.split('.')
            tmpDct = dct
            for val in splitParam[:-1]:
                tmpDct = tmpDct.setdefault(val, {})
            tmpDct[splitParam[-1]] = getattr(args, param)
    return dct


def getParamFromFullName(allParams, fullName):
    tmpDct = allParams
    for name in fullName.split('.'):
        if name not in tmpDct:
            raise KeyError('Cannot find {fullName} in existing parameters.')
        tmpDct = tmpDct[name]
    return tmpDct


def hierarchicalDictUpdate(dct1, dct2):
    for name, val2 in dct2.items():
        val1 = dct1.get(name, None)
        if val1 is None or not isinstance(val1, dict) or not isinstance(val2, dict):
            dct1[name] = val2
        else:
            hierarchicalDictUpdate(val1, val2)


def addUnknownHierarchicalArgs(allParams, addedParams, params, unknown_args):
    reParams = []
    for param in addedParams:
        parts = param.split('.')
        reStr = ''.join(
            ('[^\.]+\.' if part.startswith('sample') else f'({part}\.)?') for part in parts[:-1]) + parts[-1]
        reParams.append((re.compile(reStr), param))

    parser = argparse.ArgumentParser()
    flatParams = []
    for name, val in zip(unknown_args[::2], unknown_args[1::2]):
        if not name.startswith('--'):
            raise Exception(f'Hierarchical parameters should start with "--", got {name} instead')
        fullName = name[2:]
        if name.endswith('.__class__'):
            # Handle the special __class__ case
            tpe = utils.classFromStr
        else:
            # Normal parameters
            matchingParams = [param for pattern, param in reParams if pattern.match(fullName)]
            if len(matchingParams) == 0:
                raise ValueError(f'Cannot parse argument "{name}"')
            allTypes = set()
            for mp in matchingParams:
                defVal, tpe = getParamFromFullName(allParams, mp)
                allTypes.add(tpe)
            if len(allTypes) > 1:
                raise Exception(
                    f'Several hierarchical parameters with different types match {name}: {allTypes}.')

        addArgumentToParser(parser, fullName, tpe)
        flatParams.append(fullName)
    parsedArgs = parser.parse_args(unknown_args)
    unknown_params = buildParameterDict(parsedArgs, flatParams)

    hierarchicalDictUpdate(params, unknown_params)


if __name__ == '__main__':
    dirPath = os.path.dirname(os.path.abspath(__file__))
    sitePackagesPath = os.path.join(dirPath, '..')
    scriptPath = os.path.join(dirPath, 'blenderscript.py')

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage='python3 -m stepsblender.load [HDFPath] [-h] [--server server_address] ...',
        description='Visualize STEPS simulations with Blender.',
        epilog="""\
Optional arguments that contain dots like e.g. --Species.sampleSpecies.obj.material.color
can be used to set the properties of specific STEPS objects in the visualization.
The "sampleSpecies" part should be replaced by the STEPS name of the species in order
to get a valid argument. For example, one can set the color of species S1 with:

    --Species.S1.obj.material.color "(1, 0, 0, 1)"

But not all parts are required, it is also possible to set the same value with:

    --S1.color "(1, 0, 0, 1)"

The class of objects used for visualization can also be changed with a similar syntax.
For example, one can visualize membrane potential by using the StateDepMeshMaterial material
class:

    --Meshes.material.__class__ objects.StateDepMeshMaterial

The class needs to be given with the module it is defined in (here it is objects). Additional
parameters for setting the voltage range and colormap are listed in objects.ShaderNodeMathRescale
and objects.ShaderNodeColorMap respectively.

A higher degree of control over how data is visualized can be achieved by writing a custom
python script instead of calling this command.

For more details, see the documentation of stepsblender.HierarchicalParamReader.
        """)

    # Ressource information
    parser.add_argument(
        'HDFPath',
        action='store',
        nargs='?',
        default=argparse.SUPPRESS,
        help='Path prefix to a STEPS HDF5 file or, if the --server option is used, this should not be provided'
    )

    # Local server options
    parser.add_argument('--serverpython',
                        action='store',
                        metavar='/path/to/python',
                        help='Python binary used to launch the data loading server',
                        default='python3')

    # Blender specific options
    parser.add_argument('--blenderPath',
                        type=str,
                        action='store',
                        metavar='/path/to/blender',
                        help='Path to the blender executable',
                        default='blender')
    parser.add_argument(
        '--blenderArgs',
        type=str,
        action='store',
        metavar='/path/to/file.blend',
        help='Path to the blender file that should be used and / or other blender arguments in double quotes',
        default=None)

    # Visualization options
    allParams = HDF5BlenderLoader.listAllParameters()
    addedParams = addVisualizationParameters(parser, allParams)

    args, unknown_args = parser.parse_known_args()

    parameters = buildParameterDict(args, addedParams)
    addUnknownHierarchicalArgs(allParams, addedParams, parameters, unknown_args)

    if 'HDFPath' in args:
        if 'server' in args:
            raise ValueError('Cannot supply the --server parameter if the path to an HDF5 file was supplied.')
        servArgs = []
        for argName in ['port', 'authkey']:
            if argName in args:
                servArgs += [f'--{argName}', str(getattr(args, argName))]
        localServer = subprocess.Popen([args.serverpython, '-m', 'stepsblender.dataloader', args.HDFPath] +
                                       servArgs)
        server = 'localhost'
    elif 'server' in args:
        localServer = None
        server = args.resource

    blenderCmd = [args.blenderPath]

    if args.blenderArgs is not None:
        argsLst = args.blenderArgs.split(' ')
    else:
        argsLst = []

    if 'render' in args:
        if '-b' not in argsLst and '--background' not in argsLst:
            argsLst.append('--background')
        # Default to EXACT intersections for rendering
        parameters.setdefault('intersectAlgo', 'EXACT')

    blenderCmd += argsLst

    # Write parameter dict to file
    paramFilePath = os.path.join(tempfile.gettempdir(), '_stepsblenderParams.dat')
    with open(paramFilePath, 'wb') as f:
        pickle.dump(parameters, f)

    # Run Blender
    subprocess.run(blenderCmd + [
        '--python-use-system-env',
        '--python',
        scriptPath,
        '--',
        paramFilePath,
        sitePackagesPath,
    ])

    # Shutdown local server
    if localServer is not None:
        try:
            localServer.wait(timeout=5)
        except subprocess.TimeoutExpired:
            localServer.terminate()
