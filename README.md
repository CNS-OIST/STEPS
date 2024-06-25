Install using Docker
====================
If you don't want to do the compilation yourself and just want to quickly try STEPS, we provide a prebuilt Docker image for you. Please check https://github.com/CNS-OIST/STEPS_Docker and follow the instructions.

Install from source code
========================

To facilitate new requirements from the parallel TetOpSplit solver, STEPS 3.0 and above uses CMake system to replace the Python distutils system in previous releases. Please follow the instructions below.

Minimum Prerequisites
---------------------
1. C++ compiler supporting c++17 (e.g. gcc 7.4, clang 6)
2. Python3 (3.8 or above, 3.9 or above if using `stepsblender` Python package)
3. NumPy (http://www.numpy.org/)
4. CMake (https://cmake.org/) (3.16.3 or above)
5. Cython (http://www.cython.org/) (0.29 or above)
6. BLAS/OpenBLAS ( http://www.openblas.net/ )
7. GSL (https://www.gnu.org/software/gsl/)
8. Eigen3 (https://eigen.tuxfamily.org/)
9. METIS (https://github.com/KarypisLab/METIS)
10. `build` Python module

See install Dependencies sections

Optional Prerequisites
----------------------
1. To use the parallel SSA solver TetOpSplit: MPI libraries (e.g. MPICH https://www.mpich.org/)
2. To use the parallel EField solver: PETSc (https://www.mcs.anl.gov/petsc/)
3. To use the distributed mesh solver with bundled omega_h: Gmsh (https://gmsh.info/)

Installation From Source code
-----------------------------
_please avoid using the "Download ZIP" feature, as submodules are currently not packed in the zip file.
This includes the master branch as well as all previous releases in https://github.com/CNS-OIST/STEPS/releases._

1. Clone the repository using `git clone` in terminal, and change to the directory.
```
git clone --recursive https://github.com/CNS-OIST/STEPS.git
cd STEPS
```
(optional): To checkout a previous release, for example release with tag `4.1.0` (Release tags can be found [here](https://github.com/CNS-OIST/STEPS/tags), type in
```
git checkout tags/4.1.0 -b steps_4.1.0
git submodule update --recursive
```

2. If not already installed, install the `build` Python modules

```
pip install --user build
```

3. Run the following commands to compile the source code and install

```
git submodule update --init --recursive
mkdir build
cd build
cmake ..
make
[sudo] make install
```

Note that, by default, STEPS will install its Python dependencies during the call to `make install`, this can be prevented by giving the `-DSTEPS_INSTALL_PYTHON_DEPS=False` option to cmake.

After installation, you can check the STEPS installation with the following commands

```
python3 -c "import steps; steps._greet()"
```

If STEPS is installed successfully, you should be able to see similar information as below

```
STochastic Engine for Pathway Simulation
Version:  5.0.1
License:  GPL3.0
Website:  steps.sourceforge.net
CXX Binding: Cython
```

You can change the installation location by changing the prefix in CMake

```
cmake -DCMAKE_INSTALL_PREFIX:PATH=/MY_INSTALL_LOCATION ..
```

MPI and PETSc libraries are automatically detected in the system. If the user wants to manually choose to build STEPS with / without them it can be set

```
cmake -DUSE_MPI=[True|False] -DUSE_PETSC=[True|False] ..
```

By default, the distributed mesh solver `DistTetOpSplit` will be built. If the user wants to manually choose to build STEPS with / without this it can be set

```
cmake -DSTEPS_USE_DIST_MESH=[True|False] ..
```

By default, the distributed mesh solver will be built with bundled omega_h. An external omega_h build may be used by setting 

```
cmake -DUSE_BUNDLE_OMEGA_H=False ..
```

Please refer to [CMAKE documentation](https://cmake.org/documentation/) for customizing your installation


Simulation with serial solvers
------------------------------
STEPS 3.0 and above contain all serial solvers in previous releases, to run STEPS in serial interactive mode, open Python and import the steps module

```python
import steps
```

Scripts of serial STEPS simulations can be executed in terminal

```
python3 sim_script.py
```

Script migration in 3.6
-----------------------

A new Python API is available but scripts that worked with STEPS 3.5 should still work without any modification in STEPS 3.6.

Detailed guides for the new API can be found in the [documentation](http://steps.sourceforge.net/manual/manual_index.html).

More details in [RELEASES](./RELEASES.md) document.

Simulation with parallel solvers
--------------------------------
At the moment STEPS does not provide the interactive interface for parallel solvers `TetVesicle`, `TetOpSplit`, `DistTetOpSplit`, thus parallel simulations need to be executed via scripts in terminal e.g. with `mpirun` command

```
mpirun -n N_PROCS python3 parallel_sim_script.py
```

N_PROCS is the number of MPI processes to be created for the parallel simulation.
For the `TetVesicle` solver N_PROCS must be a minimum of 2.
For the `DistTetOpSplit` solver N_PROCS must be a power of 2.

High performance computing clusters often provide optimized parallel job scheduler and associated commands. For example, the above  `mpirun` command may be replaced by

```
srun --mpi=pmix python3 parallel_sim_script.py
```

with [Slurm](https://slurm.schedmd.com/quickstart.html) scheduler and [PMIx](https://pmix.github.io/) interface. Please refer to the documentation of your MPI solution for further customization.

Dependencies and build instructions
-----------------------------------

### Python dependencies

Since STEPS 5.0, the Python dependencies (except `build`) are automatically installed during STEPS installation. This does not apply to non-python dependencies like MPI or PETSC, which still need to be installed by the users.

On some systems, the `python3-venv` package might be needed. If so, an error message during installation will point to the exact command needed to install it on your system.

### Linux Debian based:
You can follow the installation procedure performed by the [Docker image recipe](https://github.com/CNS-OIST/STEPS_Docker/blob/329b5aaeb4f714f1bfd6fde045addbf3a338dd68/recipe/Dockerfile)

### OSX:
Use Anaconda or Miniconda:
`conda install scipy numpy matplotlib cmake cython openblas openmpi llvm-openmp`


#### Full example to build STEPS on Apple M1 Silicon with Python 3.8 through Miniconda

Prerequisites:
* Install latest XCode
* Install [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) _macOS Apple M1 64-bit_

```bash
CONDA_DIR=/path/to/miniconda
export PATH=$CONDA_DIR/bin:$PATH
conda install -c conda-forge \
  boost-cpp    \
  cmake        \
  cython       \
  eigen        \
  gfortran     \
  gmsh         \
  gsl          \
  llvm-openmp  \
  matplotlib   \
  metis        \
  mpi4py       \
  mpich        \
  numpy        \
  openblas     \
  pkg-config   \
  python-build \
  scipy
export MPICH_FC=$CONDA_DIR/bin/gfortran
export MPICH_CC=/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc
export MPICH_CXX=/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++

# Fetch and build PETSc
cd /path/to/src
git clone -b release https://gitlab.com/petsc/petsc.git petsc
pushd petsc
./configure --prefix=/path/to/petsc/installation --with-64-bit-indices=1 --with-debugging=0 \
            --with-scalar-type=real --with-x=0 --CC=mpicc --CXX=mpicxx --F77=mpif77 --FC=mpif90 \
            "MAKEFLAGS=$MAKEFLAGS"
make
make install
export PKG_CONFIG_PATH="/path/to/petsc/installation/lib/pkgconfig:$PKG_CONFIG_PATH"
popd

cd /path/to/src
git clone -b 5.0.0_beta --recursive https://github.com/CNS-OIST/STEPS.git
cd STEPS
mkdir __build
cd __build
export CC=mpicc
export CXX=mpicxx
cmake ..
make
make install
```

Validation and Examples
-----------------------
 - Short validation tests (running in a few minutes, using serial solvers only) can be found in this repository, under `test/validation`
 - Longer validation tests (using serial and parallel solvers) can be found in the [STEPS_Validation](https://github.com/CNS-OIST/STEPS_Validation) repository
 - Examples scripts (including tutorials and papers models) can be found in the [STEPS_Example](https://github.com/CNS-OIST/STEPS_Example) repository

Code Formatting and Static Analysis
-----------------------
The [hpc-coding-conventions](https://github.com/BlueBrain/hpc-coding-conventions) submodule is responsible for the code formatting and static analysis.

1. To activate code formatting of both C/C++ and CMake files, enable the CMake variable `STEPS_FORMATTING` (`-DSTEPS_FORMATTING:BOOL=ON`). This will add the following make targets: `clang-format, check-clang-format, cmake-format, check-cmake-format`.
2. To activate static analysis of C++ files with clang-tidy, enable the CMake variable `STEPS_STATIC_ANALYSIS` (`-DSTEPS_STATIC_ANALYSIS:BOOL=ON`). This will add the following make target: `clang-tidy`.

Thorough [instructions](https://github.com/BlueBrain/hpc-coding-conventions/blob/master/cpp/README.md) on how to perform code formatting and static analysis can be found in the submodule's repository.

Documentation
-------------
You can find STEPS user manual and other documentation from the STEPS official website [http://steps.sourceforge.net](http://steps.sourceforge.net)
