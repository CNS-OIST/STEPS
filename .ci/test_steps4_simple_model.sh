#!/bin/bash

set -e

# This script launches a couple of STEPS4 models
# It assumes that STEPS Python module is available in PYTHONPATH

which python
python -c "import steps"

# SSH denies key $BBPCIHPCDEPLOY_GERRIT_PRIVATE_KEY because file too permissive
install -m 400 "$BBPCIHPCDEPLOY_GERRIT_PRIVATE_KEY" "$PWD/.hpc-key"
export GIT_SSH_COMMAND="ssh -o ProxyCommand='/usr/bin/socat - PROXY:bbpproxy.epfl.ch:%h:%p' -i $PWD/.hpc-key -o IdentitiesOnly=yes"

if ! [ -d STEPS4Models ] ; then
  git clone --recursive -b main --single-branch git@github.com:CNS-OIST/STEPS4Models.git
fi
pushd STEPS4Models
git pull -r origin main

module load unstable
module load py-numpy

pushd SimpleModel/Sample_STEPS4
srun -n4 python parallel_simulation.py
popd

popd
