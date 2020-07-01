#!/bin/bash
#SBATCH --job-name=LES_
#SBATCH --partition=compute2
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --output=PF_out.log.o%j
#SBATCH --error=PF_error.log.o%j
#SBATCH --account=bb1097
#SBATCH --mail-type=FAIL


set -x


# no. of processors ntasks must be nprocx * nprocy
#ntasks=128
nprocx=256
nprocy=1
nthreads=1


# OpenMP settings
export OMP_NUM_THREADS=${nthreads}
#export MV2_ENABLE_AFFINITY=0

# limit stacksize ... adjust to your programs need
ulimit -s 102400
ulimit -c 0

#module load intel/18.0.4 
#module load openmpi/2.0.2p2_hpcx-intel14

# Settings for OpenMPI and MXM (MellanoX Messaging)
# library
export OMPI_MCA_pml=cm
export OMPI_MCA_mtl=mxm
export OMPI_MCA_mtl_mxm_np=0
export MXM_RDMA_PORTS=mlx5_0:1
export MXM_LOG_LEVEL=ERROR
# Disable GHC algorithm for collective communication
export OMPI_MCA_coll=^ghc




d_home=/pf/b/b380174/
#d_scratch=/scratch/atmomodel/kelemen/

d_nam=${d_home}PincFloit/run_Fanni/pinc_git_kelemen
fp_exe=${d_home}PincFloit/pinc/bin/pinc_master
d_work=/work/bb1097/b380174/pinc_git_kelemen/LES_imB_expl_u75_48h_h500_test4Gergo/

mkdir ${d_work}

#- go to work

cd ${d_work}

# copy namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${d_nam}/input_LES_imB_expl_u75_48h_h500_test4Gergo.f90 > input.f90

#- run the raytracer
#mpirun -np ${ntasks} ${fp_exe} 1>run.log 2>&1


srun -l --propagate=STACK --cpu_bind=cores --distribution=block:cyclic ${fp_exe} 1>run.log 2>&1

exit 0
