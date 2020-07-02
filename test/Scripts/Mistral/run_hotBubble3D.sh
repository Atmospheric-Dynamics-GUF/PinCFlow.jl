#!/bin/bash
#SBATCH --job-name=hotBubble3D
#SBATCH --partition=compute2
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:20:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097
#SBATCH --output=hotBubble3D.o%j
#SBATCH --error=hotBubble3D.e%j

# ---------------------------------------------
# Mistral batch script for the PincFloit model
# ---------------------------------------------

set -x

# no. of processors ntasks must be nprocx * nprocy
nprocx=16
nprocy=16
nproc=$((nprocx * nprocy)) 
nthreads=1
nompt=1
expname=hotBubble3D

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


d_pinc_home=/pf/b/b380792/PincFloitMSGWAM
d_nam=${d_pinc_home}/nam
fp_exe=${d_pinc_home}/src/Repo/pinc/bin/pinc
d_work=~/scratchdir/PincFloitMSGWAM/output/${expname}

mkdir -p ${d_work}

#- go to work
cd ${d_work}

# copy namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${d_nam}/NAM_${expname} > input.f90

#- run the model
srun -l --propagate=STACK,CORE --cpu_bind=cores --cpus-per-task=${nompt} ${fp_exe} 1>${d_work}/run.log 2>&1

exit 0
