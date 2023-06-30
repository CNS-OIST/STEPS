#!/bin/bash

set -e

module load unstable git python-dev

# This script launches a couple of STEPS4 models
# It assumes that STEPS Python module is available in PYTHONPATH

echo "*******************************************************************************"
echo "SPACK_INSTALLED_HASH (steps) : " ${SPACK_INSTALLED_HASH}
echo "*******************************************************************************"
which python
python -c "import steps"
echo "*******************************************************************************"

if ! [ -d STEPS4ModelRelease ] ; then
  git clone --quiet --recursive -b ${STEPS4ModelRelease_BRANCH:-main} --single-branch https://github.com/CNS-OIST/STEPS4ModelRelease.git
fi

pushd STEPS4ModelRelease
pushd SimpleModel/profiling/Sample_STEPS4/

python -m venv ./python-venv
source ./python-venv/bin/activate
pip install psutil

srun -n 4 python parallel_simulation.py

popd
popd
