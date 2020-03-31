# Distributed splitting operator

## Development environment with Spack on BB5

1. Follow installation instructions mentioned in
   [Using Spack on BB5](https://github.com/BlueBrain/spack#building-software-with-spack-at-bluebrain)
   wiki page. If Spack is already setup on your account, please ensure that
   you have the latest changes of *develop* branch.
2. Update `MODULEPATH` environment variables to take benefits of modules used by spack packages
   ```
   . /gpfs/bbp.cscs.ch/apps/hpc/jenkins/config/modules.sh
   ```
3. Clone this repository:
    ```
    git clone --recursive https://github.com/CNS-OIST/STEPS.git /path/to/steps
    pushd /path/to/steps
    git checkout distribute_steps
    ```
4. Configure spack env:
    ```
    spack setup --shebang env zee@develop
    mkdir _build
    pushd _build
    spack build-env zee@develop bash --norc
    ```
5. Build it:
    ```
    ../spconfig.py -DZee_FORMATTING:BOOL=True -DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE ..
    make -j
    ```

## Build environment on sango

```
git clone git@github.com:BlueBrain/zee.git
git checkout THE_BRANCH_TO_WORK_ON
git submodule update -r --init

STEPS_ENV=/apps/unit/DeSchutterU/steps_env/2019_06
module load $STEPS_ENV/profiles/default
module load gcc/6.2.0
module load binutils/2.27
export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$STEPS_ENV/lib/cmake"
```

Add the following line in the top CMakeLists.txt:
```diff
--- CMakeLists.txt.orig 2019-10-18 19:20:57.672038359 +0900
+++ CMakeLists.txt      2019-10-18 19:21:14.138275808 +0900
@@ -65,6 +65,7 @@
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-attributes")
   endif()
 endif()
+add_definitions("-D_GLIBCXX_USE_CXX11_ABI=0")
 message("-- CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

 if(Zee_BUILD_TESTING)
```

```
mkdir _build
cd _build
cmake \
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=True \
    -DCMAKE_C_COMPILER=icc \
    -DCMAKE_CXX_COMPILER=icpc \
    -DCMAKE_PREFIX_PATH=$STEPS_ENV \
    -DZee_OPSPLIT_CLIENTS_ONLY:BOOL=True \
    -DZee_BUILD_TESTING:BOOL=False \
    ..
make
```

## Spack integration with CLion

1. Create file `spconfig.py` with command: `spack setup --shebang env zee@develop`
1. Create file `spack-build.env` with command: `spack build-env --clean zee@develop bash --norc -c env > spack-build.env` 
1. Open zee project with CLion
1. Go to "File > Settings > Build, Execution, Deployment > ToolChains"
  * Create new "Spack" profile
  * Set "CMake" variable to the path to `spconfig.py` file
1. Go to "File > Settings > Build, Execution, Deployment > CMake"
  * Copy content of `spack-build.env` in clipboard
  * Open "Environment > ...", click on the paste button, then the "OK" one.

## Usage

### OpSplitOmega_h

* How to pass the mesh file in CLI:
  * OpSlitOmega_h: last argument of the command
* CLI Options:
  * `--test <int>`: the scenario to evaluate:
      * 0: validation test
      * 1: simple container
      * 2: multiple compartment
      * 3: diffusion only
      * 4: surface reaction validation
      * 5: a patch between two tets
      * -1 (default): simple reaction A + B => C
      * any other negative value: compare RSSA and SSA simulations of the default scenario
  * `--threshold <int>` (default 10): the diffusion operator molecules threshold
  * `--scale <float>` (default 1e-6): the import scale from the Gmsh unit to STEPS unit
  * `--num-mols-factor <int>` (default 1): multiplication factor of the initial number of molecules
  * `--debug-diffusion`: flag to track the number molecules moving between ranks
  * `--log-mesh-report`: log the mesh report before starting the simulation (default is false)
  * `--log-statedef-report`: log the statedef report before starting the simulation (default is false)
  * `--log-state-report`: log the state report before starting the simulation (default is false)
  * `--rng-seed <int>` (default is random value): force the seed of the random number generator.
    Note that previous version of OpSplitOmega_h was using the Mersenne Twister 19937 default seed value
    5489.
  * `--expected-num_reactions <int>`: Expected number of reactions during the simulation
   the program will fail if number of reactions is not in the following interval
    ```
    [
      expected_num_reactions * (1 - expected_num_reactions_tolerance),
      expected_num_reactions * (1 + expected_num_reactions_tolerance)
    ]
    ```
    Default value is -1, meaning that the number of reactions is not checked.
  * `--expected-num-reactions-tolerance <float>`: tolerance between 0 and 1
    Default value if 0.01.
  * `--expected-num-diffusions`: Expected number of reactions during the simulation
    the program will fail if number of reactions is not in
    the following interval
    ```
    [
      expected_num_reactions * (1 - expected_num_reactions_tolerance),
      expected_num_reactions * (1 + expected_num_reactions_tolerance)
    ]
    ```
    Default value is -1, meaning that the number of diffusions is not checked.
  * `--expected-num-diffusions-tolerance <float>`: Tolerance between 0 and 1
    Default value is 0.02.
  * `--do-interval <float>` (default 1e-7): step interval to use for the "diffusionOnly" scenario
  * `--end-time <float>` (default 20): simulation end time
  * `--do-inject-in-elt <int>` (default -1): For the "diffusionOnly" scenario. Tells in which element to
    inject molecules. Default is -1 meaning that molecules are spread in all elements.

## TIPS AND TRICKS

* Use `msh2osh input.msh output.osh` to convert a mesh generated with gmsh into Omega_h's owned file format
* Use the following CMake variables to control the build:
  * `Zee_CXX_OPTIMIZE`: Compile C++ with optimization, default is ON
  * `Zee_CXX_SYMBOLS`: Compile C++ with debug symbols, default is ON
  * `Zee_CXX_WARNINGS`: Compile C++ with warnings, default is ON
  * `Zee_OPSPLIT_CLIENTS_ONLY:BOOL`: Only build `OpSplitDMPlex` and `OpSplitOmega_h` clients,
    default is OFF
  * `Zee_BUILD_TESTING:BOOL`: build the tests and execute them in `make test`, default is ON

## AUTHORS

* Omar Awile <omar.awile@epfl.ch>
* Samuel Melchior <samuel.melchior@epfl.ch>
* Tristan Carel <tristan.carel@epfl.ch>
* Weiliang Chen <w.chen@oist.jp>
