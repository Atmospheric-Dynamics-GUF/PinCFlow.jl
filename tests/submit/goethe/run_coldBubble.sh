#!/bin/bash
#SBATCH --partition=general1
#SBATCH --job-name=coldBubble
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2600
#SBATCH --mail-type=FAIL
#SBATCH --time=03:00:00

set -x

# no. of processors ntasks must be nprocx * nprocy
ntasks=128
nprocx=128
nprocy=1

# OpenMP settings
export OMP_NUM_THREADS=1

# gprof for MPI
export GMON_OUT_PREFIX=gmon.out-

# env variable export
export LD_LIBRARY_PATH=LD_LIBRARY_PATH:/home/atmodynamics/jochum/spack/opt/spack/linux-scientific7-x86_64/intel-18.0.3/hypre-2.15.1-j7lo2mzfhd7bnrcivaf6bsaqjwimsp3m/lib/

dirHome=/home/atmodynamics/jochum/master
dirScratch=/scratch/atmodynamics/jochum/master

dirNam=${dirHome}/pincflow_tfc_semi_implicit_boussinesq/tests/input
exe=${dirHome}/pincflow_tfc_semi_implicit_boussinesq/bin/pinc
dirWork=${dirScratch}/pincflow_tfc_semi_implicit_boussinesq/tests/coldBubble

mkdir ${dirWork}

# go to work

cd ${dirWork}
rm *

# copy namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${dirNam}/input_coldBubble.f90 > input.f90

# run the raytracer
mpirun -np ${ntasks} ${exe} 1>run.log 2>&1

exit 0
