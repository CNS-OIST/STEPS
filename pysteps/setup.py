#!/usr/bin/env python
from __future__ import print_function

"""A setuptools based setup module."""
from setuptools import setup, find_packages, Extension
import os, sys, sysconfig
from os import path
from distutils.sysconfig import get_config_var
from Cython.Build import cythonize

cur_path = path.abspath(path.dirname(__file__))
cflags = get_config_var('CFLAGS')
os.environ['CPPFLAGS'] = cflags + ' ' + '-std=c++11'
lib_dirs = os.environ.get('LD_LIBRARY_PATH','').split(':')


def distutils_dir_name():
    """Returns the name of a distutils build directory"""
    f = "lib.{platform}-{version[0]}.{version[1]}"
    return f.format(platform=sysconfig.get_platform(),
                    version=sys.version_info)

# ----------------------------------------------
# Command line args ##
# ----------------------------------------------
# Just get build-dir
if '--print-build-dir' in sys.argv:
    print(distutils_dir_name())
    exit(0)

# LD Libraries - No spaces please
for opt in sys.argv[:]:
    if opt.startswith('-L'):
        lib_dirs.append(opt[2:])
        sys.argv.remove(opt)


# With mpi?
NO_MPI = '--nompi' in sys.argv
if NO_MPI:
    sys.argv.remove('--nompi')


# ----------------------------------------------
# Extension setup
# ----------------------------------------------
steps_ext = Extension('steps.cysteps_mpi', ['pysteps/cysteps_mpi.pyx'], include_dirs = [ 'src' ], libraries=['steps', 'mpich'], library_dirs=lib_dirs, language='c++') if not NO_MPI \
       else Extension('steps.cysteps',     ['pysteps/cysteps.pyx'],     include_dirs = [ 'src' ], libraries=['steps'],          library_dirs=lib_dirs, language='c++')

extensions = [
    steps_ext
  ]

setup_opts = {
    'name'         : 'steps',
    'version'      : '3.1.0',
    'author'       : 'STEPS Development Team',
    'author_email' : 'steps.dev@gmail.com',
    'url'          : 'http://steps.sourceforge.net',
    'description'  : 'STochastic Engine for Pathway Simulation',
    'download_url' : 'http://sourceforge.net/projects/steps/files/src',
    'platforms'    : ['Mac OS X', 'Windows XP', 'Windows Vista', 'Linux', 'Unix'],
    'license'      : 'GNU General Public License Version 3.0',

    'package_dir'  : {'':'pysteps'},
    'packages'     : find_packages('pysteps'),
    'package_data' : {},
    'ext_modules'  : cythonize(extensions, include_path=['.'], gdb_debug=False),
}

if __name__ == '__main__':
    print('Lib Path:', lib_dirs)
    setup( **setup_opts )
