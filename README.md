Install using Docker
====================
If you don't want to do the compilation yourself and just want to quickly try STEPS, we provide a prebuilt Docker image for you. Please check https://github.com/CNS-OIST/STEPS_Docker and follow the instructions.

Install from source code
========================

To facilitate new requirements from the parallel TetOpSplit solver,
STEPS 3.0 and above uses CMake system to replace the Python distutils system
in previous releases. Please follow the instructions below.

Minimum Prerequisites
---------------------
1. C++ compiler supporting c++11 (e.g. gcc 4.8, clang 3.3)
2. Python2/3 (2.7.x / 3.3.x or above)
3. NumPy (http://www.numpy.org/)
4. CMake (https://cmake.org/)
5. Cython (http://www.cython.org/)
6. BLAS/OpenBLAS ( http://www.openblas.net/ )

See install Dependencies sections

Optional Prerequisites
----------------------
1. To use the parallel SSA solver TetOpSplit: MPI libraries (e.g. MPICH https://www.mpich.org/)
2. To use the parallel EField solver: PETSc (https://www.mcs.anl.gov/petsc/)


Removing Previous STEPS Installation
------------------------------------
Due to potential conflict, please remove any previous installed version
below version 3.0.0, to check and remove previous release,
enter the following commands in a terminal:

```
python -c "import steps; print steps.__file__"
```

The commands will provide you the path of your current STEPS installation, for example

```
/usr/local/lib/python2.7/site-packages/steps/__init__.pyc
```

you can remove the STEPS installation in the above location via terminal command:

```
[sudo] rm -rf /usr/local/lib/python2.7/site-packages/steps
```

Installation From Source code
-----------------------------
_please avoid using the "Download ZIP" feature, as submodules are currently not packed in the zip file.
This includes the master branch as well as all previous releases in https://github.com/CNS-OIST/STEPS/releases._

1. Clone the repository using `git clone` in terminal, and change to the directory.
```
git clone --recursive https://github.com/CNS-OIST/STEPS.git
cd STEPS
```
(optional): To checkout a previous release, for example release with tag `3.4.0` (Release tags can be found [here](https://github.com/CNS-OIST/STEPS/tags), type in
```
git checkout tags/3.4.0 -b steps_3.4.0
git submodule update --recursive
```
2. run the following commands to compile the source code and install

```
git submodule update --init --recursive
mkdir build
cd build
cmake ..
make
[sudo] make install
```

You can change the installation location by changing the prefix in CMake

```
cmake -DCMAKE_INSTALL_PREFIX:PATH=/MY_INSTALL_LOCATION ..
```

MPI and PETSc libraries are automatically detected in the system. If the user
wants to manually choose to build STEPS with / without them it can set

```
cmake -DUSE_MPI=[True|False] -DUSE_PETSC=[True|False] ..
```

Please refer to [CMAKE documentation](https://cmake.org/documentation/) for customizing your installation


After installation, you can check the STEPS installation with the following commands

```
python -c "import steps; steps._greet()"
```

If STEPS is installed successfully, you should be able to see similar information as below

```
STochastic Engine for Pathway Simulation
Version:  3.4.0
License:  GPL2.0
Website:  steps.sourceforge.net
CXX Binding: Cython
```

Simulation with serial solvers
------------------------------
STEPS 3.0 and above contain all serial solvers in previous releases,
to run STEPS in serial interactive mode, open Python and import the steps module

```python
import steps
```

Scripts of serial STEPS simulations can be executed in terminal

```
python sim_script.py
```

Script migration in 3.5
-----------------------

The data type of mesh entities like tetrahedrons and triangles has changed in STEPS 3.5. There are
now unsigned integers and a new CMake option controls their size: 32 bits (the default) or 64 bits.
A few macros have been introduced in the Python API to keep the simulation code independent of
the STEPS installation.

1. Use `steps.geom.UNKNOWN_TET` and `steps.geom.UNKNOWN_TRI` instead of -1, for instance:

   ```python
   if mesh.findTetByPoint(a_point) == steps.UNKNOWN_TET:
       raise RuntimeError("Point is not in the mesh")
   ```

   instead of

   ```python
   if mesh.findTetByPoint(a_point) == -1:
       raise RuntimeError("Point is not in the mesh")
   ```

A Flake8 extension allows you to track these cases in your code, see [RELEASES]((./RELEASES.md)) document.

1. Use `steps.geom.INDEX_DTYPE` to instantiate NumPy arrays. for instance:

   ```python
   tetrahedrons = np.array(num_tets, dtype=steps.geom.INDEX_DTYPE)
   ```

   instead of

   ```python
   tetrahedrons = np.array(num_tets, dtype=np.uint32)
   ```

More details in [RELEASES](./RELEASES.md) document.

Simulation with parallel TetOpSplit
-----------------------------------
At the moment STEPS does not provide the interactive interface for parallel TetOpSplit solver,
thus parallel simulations need to be executed via scripts in terminal with "mpirun" command

```
mpirun -n N_PROCS python parallel_sim_script.py
```

N_PROCS is the number of MPI processes to be created for the parallel simulation.

Please refer to the documentation of your MPI solution for further customization.


Dependencies
-------------
Linux Debian based:
 `apt-get install g++ gcc cmake libopenblas-dev libmpich-dev python-numpy cython python-scipy`

OSX:
 `brew install openblas mpich`

If you use Anaconda:
`conda install scipy numpy matplotlib cmake cython openblas openmpi`

Validation and Examples
-----------------------
 - Short validation tests (running in a few minutes, using serial solvers only) can be found in this repository, under `test/validation`
 - Longer validation tests (using serial and parallel solvers) can be found in the [STEPS_Validation](https://github.com/CNS-OIST/STEPS_Validation) repository
 - Examples scripts (including tutorials and papers models) can be found in the [STEPS_Example](https://github.com/CNS-OIST/STEPS_Example) repository

Documentation
-------------
You can find STEPS user manual and other documentation from the STEPS official website [http://steps.sourceforge.net](http://steps.sourceforge.net)
