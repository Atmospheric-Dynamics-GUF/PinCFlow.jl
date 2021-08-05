#!/bin/bash
##general1/general2 = new/old cores
#SBATCH --partition=test

#SBATCH --job-name=GW2sim2

##not udes anymore ...
##SBATCH --constraint=intel20

##switch off multi-threading
##BATCH --hint=nomultithread

#SBATCH --ntasks=1

##usually nodes determined by machine. However, if more memory is needed ...
##SBATCH --nodes=2

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2600
#SBATCH --mail-type=FAIL
#SBATCH --time=00:20:00

set -x

#echo "------------------------------------------------------------"
#echo LOADL_STEP_ID=$LOADL_STEP_ID
#echo "------------------------------------------------------------"

#set -x

# no. of processors ntasks must be nprocx * nprocy
ntasks=1
nprocx=1
nprocy=1

# OpenMP settings
export OMP_NUM_THREADS=1

# gprof for MPI
export GMON_OUT_PREFIX=gmon.out-

# env variable export
export LD_LIBRARY_PATH=LD_LIBRARY_PATH:/home/atmodynamics/schmid/spack/opt/spack/linux-scientific7-x86_64/intel-18.0.3/hypre-2.15.1-j7lo2mzfhd7bnrcivaf6bsaqjwimsp3m/lib/


d_home=/home/atmodynamics/schmid/git_PF_2020/pinc_f2Omegasiny_BK19cases
d_scratch=/scratch/atmodynamics/schmid/paper2020_results
#d_home2=/home/atmodynamics/schmid/git_PF_2020/pinc
#d_home2=/home/atmodynamics/schmid/test_PFApril2020

d_nam=${d_home}/bin_IGW2_si_muscl2
fp_exe=${d_home}/bin_IGW2_si_muscl2/pinc
d_work=${d_scratch}/BK19_IGWcase2_si_muscl2

mkdir ${d_work}

#- go to work

cd ${d_work}

# copy namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${d_nam}/input.f90 > input.f90

#- run the raytracer
mpirun -np ${ntasks} ${fp_exe} 1>run.log 2>&1

exit 0
