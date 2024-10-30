#!/bin/env bash

#SBATCH --job-name=whatevs
#SBATCH --account=nn9849k
#SBATCH --time=1-0:0:0
#SBATCH --ntasks=1 --cpus-per-task=2 --mem-per-cpu=8GB
#SBATCH --array=1-10

set -o errexit # exit on errors
set -o nounset # treat unset variables as errors

lakes="Ieŝjávri foobar whatevs"

lake=$(echo $lakes|cut -d' ' -f$SLURM_ARRAY_TASK_ID)

for run in {1..20}; do
    ./julia.sif demo.jl $lake $run
done