#!/bin/env bash

#SBATCH --job-name=something-something-identifying
#SBATCH --account=nn9849k
#SBATCH --time=1-0:0:0
#SBATCH --ntasks=1 --cpus-per-task=2 --mem-per-cpu=8GB
#SBATCH --array=1-20

### USAGE, e.g.:
###
### for lake in Ieŝjávri lake2 lake3; do
###     sbatch --job-name=$lake ./batchsubmit.sh $lake
### done
###
###

set -o errexit # exit on errors
set -o nounset # treat unset variables as errors

lake=$1 # name of lake passed as first argument

./julia.sif demo.jl $lake $SLURM_ARRAY_TASK_ID
