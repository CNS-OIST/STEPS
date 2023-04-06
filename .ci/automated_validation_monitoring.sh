#!/bin/bash

set -e

module load unstable git python-dev py-notebook

python -m venv venv
source ./venv/bin/activate
pip install seaborn
pip install pyarrow

echo "********************************************************************"
echo "   Useful locations"
echo "********************************************************************"

base_dir=/gpfs/bbp.cscs.ch/data/project/proj12/jenkins/subcellular/HBP_STEPS/automated_validation_monitoring
date=$(git show -s --format=%cd --date=short)
shorthash=$(git show -s --format=%h)
export res_dir=${base_dir}/${model_name}/PR_${CI_EXTERNAL_PULL_REQUEST_IID}_${date}_STEPS${steps_version}_${shorthash}
export ref_dir=${base_dir}/${model_name}/${compare_with}

echo base_dir:
echo ${base_dir}
echo res_dir:
echo ${res_dir}
echo ref_dir:
echo ${ref_dir}
cd ..
echo curr_dir:
echo $(pwd)

echo "********************************************************************"
echo "   Test python and steps"
echo "********************************************************************"

which python
python -c "import steps"

# SSH denies key $BBPCIHPCDEPLOY_GERRIT_PRIVATE_KEY because file too permissive
install -m 400 "$BBPCIHPCDEPLOY_GERRIT_PRIVATE_KEY" "$PWD/.hpc-key"
export GIT_SSH_COMMAND="ssh -i $PWD/.hpc-key -o IdentitiesOnly=yes"

echo "********************************************************************"
echo "   Clone STEPS4ModelRelease"
echo "********************************************************************"

git clone --recursive -b ${STEPS4ModelRelease_BRANCH:-main} --single-branch https://github.com/CNS-OIST/STEPS4ModelRelease.git

pushd STEPS4ModelRelease/${simulation_dir}

echo "********************************************************************"
echo "   Run ${model_name} simulations"
echo "********************************************************************"

ACCOUNT_=`sacct --format=account -j $SLURM_JOBID | tail -n1`

echo "ACCOUNT_" ${ACCOUNT_//[[:blank:]]/}
echo "SLURM_JOB_USER" $SLURM_JOB_USER
echo "SLURM_JOB_PARTITION" $SLURM_JOB_PARTITION

params="--account=${ACCOUNT_//[[:blank:]]/} --partition=${SLURM_JOB_PARTITION} --constraint=cpu run.sbatch"

echo "PARAMS: " ${params}

echo "********************************************************************"
echo "   Submit simulations"
echo "********************************************************************"

jobids=()
jobids=(${jobids[@]} $(sbatch --parsable --job-name=${slurm_job_name} $params))

echo "********************************************************************"
echo "   Jobs submitted"
echo "   Waiting for the jobs"
echo "********************************************************************"

jobs_pending=1
while [[ $jobs_pending != 0 ]]
do
    sleep 60

    jobs_pending=$(squeue -h -j$(echo ${jobids[@]} | sed -z 's/ /,/g') 2>/dev/null | wc -l)
    echo "jobs_pending=$jobs_pending"
done

echo "********************************************************************"
echo "   Jobs completed"
echo "   Clone STEPS_Validation"
echo "********************************************************************"

# if STEPS_Validation is missing clone only the postproc folder
if ! [ -d STEPS_Validation ] ; then
  mkdir STEPS_Validation
  pushd STEPS_Validation
  git init
  git remote add -f origin https://github.com/CNS-OIST/STEPS_Validation.git
  git config core.sparseCheckout true
  echo 'postproc' >> .git/info/sparse-checkout
  git pull origin master
  popd
fi

echo "********************************************************************"
echo "   Get the notebook"
echo "********************************************************************"

# get the nb. We will need it later
popd
mv STEPS4ModelRelease/validation/automated_validation_monitoring_results_displayer.ipynb STEPS4ModelRelease/${simulation_dir}/STEPS_Validation/postproc
#go in the folder to run the nb
pushd STEPS4ModelRelease/${simulation_dir}/STEPS_Validation/postproc

echo "********************************************************************"
echo "   Run postproc as in the example ${model_name}.py and pretty record results "
echo "********************************************************************"

## execute the jupyter notebook
jupyter-nbconvert --execute --to html --no-input automated_validation_monitoring_results_displayer.ipynb

echo "********************************************************************"
echo "   Print html "
echo "********************************************************************"

cat automated_validation_monitoring_results_displayer.html

echo "********************************************************************"
echo "   Copy refined data and pictures to: ${res_dir} "
echo "********************************************************************"

mkdir -p ${res_dir}/raw_traces/STEPS${steps_version}/
setfacl -R -m g:bbp-user-proj16_AppGrpU:wX ${res_dir}
cp -r ../../raw_traces/STEPS${steps_version}/res* ${res_dir}/raw_traces/STEPS${steps_version}/
setfacl -R -m g:bbp-user-proj16_AppGrpU:wX ${res_dir}
cp -r ${model_name}/pics ${res_dir}/
setfacl -R -m g:bbp-user-proj16_AppGrpU:wX ${res_dir}
cp automated_validation_monitoring_results_displayer.html ${res_dir}
setfacl -R -m g:bbp-user-proj16_AppGrpU:wX ${res_dir}

