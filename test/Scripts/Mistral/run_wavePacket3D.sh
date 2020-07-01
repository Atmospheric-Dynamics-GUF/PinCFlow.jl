#!/bin/bash

#SBATCH --job-name=testcase1
#SBATCH --partition=compute2
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:30:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097
#SBATCH --output=testcase1.o%j
#SBATCH --error=testcase1.e%j


# ---------------------------
# vTool script for Goethe HLR
# ---------------------------

# cpu settings
nprocx=128
nprocy=4
nproc=$((nprocx * nprocy))
nthreads=1
nompt=1

# load modules
module load intel/17.0.6
module load openmpi/2.0.2p2_hpcx-intel14

# Bind your OpenMP threads
export OMP_NUM_THREADS=${nthreads}
export KMP_AFFINITY=verbose,granularity=core,compact,1
export KMP_STACKSIZE=64m

# settings for OpenMPI and MXM (MellanoX Messaging) library
export OMPI_MCA_pml=cm
export OMPI_MCA_mtl=mxm
export OMPI_MCA_mtl_mxm_np=0
export MXM_RDMA_PORTS=mlx5_0:1
export MXM_LOG_LEVEL=ERROR
# disable GHC algorithm for collective communication
export OMPI_MCA_coll=^ghc

# limit stacksize ... adjust to your programs need
ulimit -s 102400
ulimit -c 0

# path for the binary
fp_exe=$PWD/../../../pinc/bin/pinc
[[ -x ${fp_exe} ]] || { echo 'Binary file does not exist!' ; exit ; }

# namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" NAM > input.f90

# run the binary
srun -l --propagate=STACK,CORE --cpu_bind=cores --cpus-per-task=${nompt} ${fp_exe} 1>run.log 2>&1

# create case summary
. create-summary
