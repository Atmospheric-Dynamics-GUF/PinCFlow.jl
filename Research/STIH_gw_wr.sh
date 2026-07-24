#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=STIHwr
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --hint=nomultithread
#SBATCH --mem=220G
#SBATCH --time=0-02:00:00
#SBATCH --mail-type=ALL
#SBATCH --account=bb1097

set -x

# Set Intel MPI configuration on compute partition.
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
srun --distribution=block:cyclic julia --project=Research -e 'include("Research/Research.jl"); using .Research; wp(; 
k = 2*pi/30.0e3, 
l = 0.0, 
m = 2*pi/3.0e3, 
branch = -1,
version = 2.0, 
rx = 0.0, 
ry = 0.0, 
rz = 5.0e3,
a0 = 0.5,
coriolis_frequency = 0.0,
lx = 30.0e3,
ly = 30.0e3,
lz = 80.0e3,
x0 = 0.0,
y0 = 0.0,
z0 = 10.0e3,
x_size = 32,
y_size = 1,
z_size = 854,
npx = 8,
output_file = "STIH_gw_wr.h5",
output_interval = 360.0, 
tmax = 3600.0 * 15, 
turbulence_scheme = :TKEScheme, 
tracer_setup = :NoTracer,
prepare_restart = false,
momentum_coupling = true,
entropy_coupling = true,
turbulence_filter_order = 0,
restart = false,
input_file = "STIH_gw_wr.h5")' &> STIH_gw_wr.log

exit 0