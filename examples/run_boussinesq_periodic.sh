#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=boussinesq_periodic
#SBATCH --nodes=14
#SBATCH --ntasks-per-node=8
#SBATCH --hint=nomultithread
#SBATCH --mem=24G
#SBATCH --time=0-08:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

# Set Intel MPI configuration on compute partition.
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
srun --distribution=block:cyclic julia --project=examples -e 'using PinCFlow; 
for a0 in 0.6:0.2:1.6
    for lambda in 50:50:300
        output_file = "boussinesq_periodic_a0-$(Int(a0*10))_lambda-$(lambda).h5"

        boussinesq_periodic(;
            npx = 28,
            npy = 4,
            npz = 1,
            output_file = output_file,
            lambdah = lambda,
            a0 = a0,
        )
    end
end' &> boussinesq_periodic.log

exit 0
