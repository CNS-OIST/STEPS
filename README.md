STEPS 3.0 Instructions
======================

To facilitate new requirements from the parallel TetOpSplit solver,
STEPS 3.0 uses CMake system to replace the Python distutil system
in previous releases. Please follow the instructions below.

Prerequisites
-------------
1. C++ compiler such as g++, if you weant to enable the parallel solver,
an MPI compiler such as OpenMPI is also required.

2. Python (We recommand Python 2.7.x)
3. NumPy (http://www.numpy.org/)
4. CMake (https://cmake.org/)
5. SWIG (http://www.swig.org/)

Removing Previous STEPS Installation
------------------------------------
Due to potential conflict, please remove any previous installed version
before the installation of STEPS 3.0, to check and remove previous release,
enter the following commands in a terminal:

        python -c "import steps; print steps.__file__"

The commands will provide you the path of your current STEPS installation, for example

       /usr/local/lib/python2.7/site-packages/steps/__init__.pyc

you can remove the STEPS installtion in the above location via terminal command:

       [sudo] rm -rf /usr/local/lib/python2.7/site-packages/steps

Installation From Source code
-----------------------------
1. Clone the repository using the folowing command in terminal

        git clone https://github.com/CNS-OIST/STEPS_Release.git

2. run the following commands to compile the source code and install

        cd STEPS_Release
        mkdir build
        cd build
        cmake ..
        make
        [sudo] make install

You can change the installation location by changing the prefix in CMake

        cmake -DCMAKE_INSTALL_PREFIX:PATH=/MY_INSTALL_LOCATION ..

Please refer to cmake documentation for customizing your installation
https://cmake.org/documentation/

Simulation with serial solvers
------------------------------
STEPS 3.0 contains all serial solvers in previous releases,
to run STEPS 3.0 in serial interactive mode, open Python and import the steps module

        import steps

Scripts of serial STEPS simulations can be exexuted in terminal

        python sim_script.py

Simulation with parallel TetOpSplit
-----------------------------------
At the moment STEPS does not provide the interactive interface for parallel TetOpSplit solver,
thus parallel simulations need to be executed via scripts in terminal with "mpirun" command

        mpirun -n N_PROCS python parallel_sim_script.py

N_PROCS is the number of MPI processes to be created for the parallel simulation.

Please refer to the documentation of your MPI solution for further customization.

