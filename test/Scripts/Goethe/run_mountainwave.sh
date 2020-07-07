#!/bin/bash
#SBATCH --partition=general1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2600
#SBATCH --mail-type=FAIL
#SBATCH --time=02:30:00
#SBATCH --extra-node-info=2:20:1

set -x

# load MPI module (openmpi)
#module load openmpi/intel-17.0.1/1.8.1

#echo "------------------------------------------------------------"
#echo LOADL_STEP_ID=$LOADL_STEP_ID
#echo "------------------------------------------------------------"

#set -x

# no. of processors ntasks must be nprocx * nprocy
ntasks=128
nprocx=128
nprocy=1

# OpenMP settings
export OMP_NUM_THREADS=1
#export MV2_ENABLE_AFFINITY=0

d_home=/home/atmomodel/kelemen/
d_scratch=/scratch/atmomodel/kelemen/

d_nam=${d_home}run_Fanni/pinc_new_rec1_topo
fp_exe=${d_home}PincFloit/PincFloit_new_rec1_topo/pinc
d_work=${d_scratch}pinc_new_rec1_topo/LES_imB_expl_u75_48h_h500/

mkdir ${d_work}

#- go to work

cd ${d_work}

# copy namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${d_nam}/input_LES_imB_expl_u75_48h_h500.f90 > input.f90

#- run the raytracer
mpirun -np ${ntasks} ${fp_exe} 1>run.log 2>&1

exit 0
