#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=wp_1d
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --mem=90G
#SBATCH --time=0-08:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

# Set Intel MPI configuration on compute partition.
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
srun --distribution=block:cyclic julia --project=Research -e 'include("Research/Research.jl"); using .Research; wp(; x_size=32, y_size=1, z_size=854, npx=4, lx=30.e3, ly=30.e3, lz=80.e3, k=2*pi/30.e3, l=0.0, m=2*pi/3.e3, rx=0., ry=0., rz=5.e3, x0=0., y0=0., z0=10.e3, a0=0.5, version=2.0, coriolis_frequency=0., output_file="boeloeni_STIH_wr.h5", output_interval=15*3600.0/100.0, tmax=15*3600.0, tracer_setup=:NoTracer)' &> boeloeni_STIH_wr.log

exit 0