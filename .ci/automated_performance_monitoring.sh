#!/bin/bash

set -e

module load unstable git python python-dev py-notebook

echo "*******************************************************************************"
echo "SPACK_INSTALLED_HASH (steps) : " ${SPACK_INSTALLED_HASH}
echo "*******************************************************************************"
which python
python -c "import steps"
echo "*******************************************************************************"

# SSH denies key $BBPCIHPCDEPLOY_GERRIT_PRIVATE_KEY because file too permissive
install -m 400 "$BBPCIHPCDEPLOY_GERRIT_PRIVATE_KEY" "$PWD/.hpc-key"
export GIT_SSH_COMMAND="ssh -i $PWD/.hpc-key -o IdentitiesOnly=yes"

###############################################################################

comparison_dir=/gpfs/bbp.cscs.ch/data/project/proj12/jenkins/subcellular/HBP_STEPS/automated_performance_monitoring
export comparison_with=$comparison_dir/$PERF_MON_MODEL_NAME/$compare_with

subcellular_meshes=/gpfs/bbp.cscs.ch/data/project/proj12/jenkins/subcellular/HBP_STEPS/meshes
if [[ ${PERF_MON_MODEL_NAME} = CaBurstBackground ]] || [[ ${PERF_MON_MODEL_NAME} = CaBurstFullModel ]]
then
  export CaBurstUnifiedMesh=$subcellular_meshes/$mesh
fi

if ! [ -d STEPS4ModelRelease ] ; then
  git clone --quiet --recursive -b ${STEPS4ModelRelease_BRANCH:-main} --single-branch https://github.com/CNS-OIST/STEPS4ModelRelease.git
fi

pushd STEPS4ModelRelease

if [[ ${PERF_MON_MODEL_NAME} = SimpleModel ]]
then
  pushd SimpleModel/profiling
elif [[ ${PERF_MON_MODEL_NAME} = CaBurstBackground ]]
then
  pushd caBurst/background/profiling
elif [[ ${PERF_MON_MODEL_NAME} = CaBurstFullModel ]]
then
  pushd caBurst/full/profiling
else
  echo "PERF_MON_MODEL_NAME is not recognized:  ${PERF_MON_MODEL_NAME}. "
  echo " Use a known model for the performance monitoring."
  exit 1
fi

# execute the jupyter notebook and save the output as html file
jupyter-nbconvert \
--execute \
--to html \
--no-input \
CI_Perf.ipynb

echo "Print html *******************************************************************************"
cat CI_Perf.html
echo "******************************************************************************* Print html"

popd
popd
