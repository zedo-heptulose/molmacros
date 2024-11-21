#!/bin/bash

#SBATCH --job-name=HBPO_tbu
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p genacc_q
#SBATCH -t 0-00:05:00
#SBATCH --mem-per-cpu=2GB

conda init bash
source ~/.bashrc
conda activate crest_3
export OMP_STACKSIZE=4G
export OMP_MAX_ACTIVE_LEVELS=1
ulimit -s unlimited
xtb HBPO_tbu.xyz > HBPO_tbu.out --gfn ff --chrg 0 --uhf 0 --opt --cma
