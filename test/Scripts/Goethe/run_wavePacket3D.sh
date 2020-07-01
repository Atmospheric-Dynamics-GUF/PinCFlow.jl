#!/bin/bash

#SBATCH --partition=general2
#SBATCH --ntasks=512
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2600
#SBATCH --mail-type=FAIL
#SBATCH --time=00:30:00

# ---------------------------
# vTool script for Goethe HLR
# ---------------------------

# cpu settings
ntasks=512
nprocx=128
nprocy=4
threads=1

# openmp settings
export OMP_NUM_THREADS=${threads}

# path for the binary
fp_exe=$PWD/../../../pinc/bin/pinc
[[ -x ${fp_exe} ]] || { echo 'Binary file does not exist!' ; exit ; }

# namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" NAM > input.f90

# run the binary
/usr/bin/time -v mpirun -np ${ntasks} ${fp_exe} 1>run.log 2>&1

# create case summary
. create-summary
