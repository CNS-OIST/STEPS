# Testing

This document contains the STEPS code testing guidelines. It should answer any
questions you may have as an aspiring STEPS contributor.

## Test suites

STEPS has four test suites:

* C++ unit tests located in `test/unit/cpp`
* Python unit tests located in `test/unit/py`
* Minimal integration tests in `test/validation`
* Validation tests located in [CNS-OIST/STEPS_Validation](https://github.com/CNS-OIST/STEPS_Validation)
repository.

## Writing new tests

Most code changes will fall into one of the following categories.

### Writing tests for new features

New code should be covered by unit tests. If the code is difficult to test with
a unit tests then that is a good sign that it should be refactored to make it
easier to reuse and maintain.

### Writing tests for bug fixes

Bugs fixes should include a unit test case which exercises the bug.

A bug fix may also include new assertions in an existing integration tests for the
API endpoint.

### Validation tests

FIXME: define guidelines

## Running tests

### Unit tests
To run the unit test suite:

```
cd /path/to/build/path
make check
```

To execute only one test: `ctest --output-on-failure -R TEST_NAME`

### Integration tests

To run the integration test suite, follow instructions of
[CNS-OIST/STEPS_Validation](https://github.com/CNS-OIST/STEPS_Validation)
repository.


## Code Analysis

### Code coverage

Code coverage provides report to know how much the STEPS code is executed
by unit tests. Analysis provides detailed reports in both HTML and XML
showing the number of times every line of code is executed when running
unit tests.

### Requirements

Following tools are required: gcov, lcov, gcovr

#### Run coverage

Compilation with gcc is required to perform code coverage.

```
cd /path/to/STEPS
mkdir _build && pushd _build
cmake -DENABLE_CODECOVERAGE:BOOL=TRUE -DUSE_MPI:BOOL=False -DUSE_PETSC:BOOL=False -DCMAKE_BUILD_TYPE=Debug ..
make -j all coverage_init test coverage
```

Under the hood:

1. `all`: compilation phase is different: every object file is attached to a .gcno file,
containing information to reconstruct the basic block graphs and assign source line
numbers to blocks.
1. `coverage_init`: creates .gcda count data file for every .gcno,
with counters reset to 0
1. `test`: counters in .gcda files are updating while executing the test suite.
1. `coverage`: generate both XML and HTML reports based on .gcno and .gcda

### Valgrind

Valgrind is an instrumentation framework to build analysis from running applications.


#### Instrument unit test suite

It is strongly advise to run C++ unit-tests with valgrind memory checker activated.

```
cd /path/to/STEPS
mkdir _build && pushd _build
cmake -DCMAKE_BUILD_TYPE=Debug -DVALGRIND=valgrind ..
export VALGRIND_OPTS="--tool=memcheck --track-origins=yes --leak-check=full --show-leak-kinds=all --verbose"
make test
```

#### Selectively suppress errors

Simply add suppression files to `test/ci/valgrind` directory with *.supp* extension
to manually remove some warnings.

#### Generate suppression files

To generate a suppression file:
1. Instrument your program with the following valgrind option:
`--gen-suppressions=all --log-file=report.log`
1. Extract all `{ }` sections from the log file your want to permanently ignore
1. Put them in *.supp* files in `test/ci/valgrind/` directory. Typically, you should
   have one file per library or namespace.

CMake won't detect that `test/ci/valgrind` has been modified, so it is required
to run `cmake .` explicitly after adding or removing a *.supp* file.

### Sanitizers

#### ASAN on Blue Brain 5

1. setup the required build environment, for instance with spack: `spack build-env steps@develop bash --norc`
2. Load the most recent LLVM module: `module load unstable llvm`
   Ensure the `CC` and `CXX` environment variables point to this compiler.
3. Configure STEPS with the sanitizer enabled:
   ```bash
   mkdir build-asan ; pushd build-asan
   cmake -DSTEPS_SANITIZERS:STRING="address" \
         -DPYTHON_EXECUTABLE=`which python3` \
         -DUSE_BUNDLE_SUNDIALS:BOOL=FALSE \
         -DUSE_BUNDLE_OMEGA_H:BOOL=FALSE ..
   ```

   * `STEPS_SANITIZERS`: comma separated list of checks, see
      https://github.com/BlueBrain/hpc-coding-conventions/blob/master/cpp/cmake/sanitizers.cmake#L1-L6
   * `PYTHON_EXECUTABLE`: needs to be forced otherwise `/usr/bin/python3` is picked
   * `USE_BUNDLE_*`: to speed-up the compilation

4. Build STEPS: `make -j`
5. Execute the tests with sanitizer enabled: `ctest`
   To execute an arbitrary Python script using STEPS, mimic the environment variables set by ctest
   to run the tests, see `ctest -VV` to have the full command written to the standard output.

Known issues:
* with srun/ASAN on BB5: `MPT ERROR: Unable to restrict trim threshold for heap.`
  The solution is to set the environment variable: `SUPPRESS_XPMEM_TRIM_THRESH=1`
