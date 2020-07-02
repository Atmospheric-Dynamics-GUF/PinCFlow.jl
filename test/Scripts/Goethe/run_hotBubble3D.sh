#!/bin/bash
#SBATCH --partition=general2
#SBATCH --ntasks=256
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-type=FAIL
#SBATCH --time=00:10:00

# ----------------------------------------------------------------------------
# Goethe Cluster batch script for PincFloitMSGWAM
# ----------------------------------------------------------------------------

set -x

# exp name
case=hotBubble3D

# cpu settings
ntasks=256
nprocx=16
nprocy=16
threads=1

# openmp settings
export OMP_NUM_THREADS=${threads}

# path for work dir
d_work=/scratch/atmodynamics/boeloeni/PincFloitMSGWAM

# path for the exp
d_exp=${d_work}/output/${case}

# path for the namelist
d_nam=${HOME}/PincFloitMSGWAM/nam

# path for the binary
fp_exe=${HOME}/PincFloitMSGWAM/src/Repo/pinc/bin/pinc

# go to exp dir
if [ ! -d ${d_exp} ]; then
    mkdir -p ${d_exp}
fi
cd ${d_exp}

# copy namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${d_nam}/NAM_${case} > input.f90

# run the binary
/usr/bin/time -v mpirun -np ${ntasks} ${fp_exe}
