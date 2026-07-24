#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=igw_wr
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=8
#SBATCH --hint=nomultithread
#SBATCH --mem=220G
#SBATCH --time=0-08:00:00
#SBATCH --mail-type=ALL
#SBATCH --account=bb1097

set -x

# Set Intel MPI configuration on compute partition.
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
srun --distribution=block:cyclic julia --project=Research -e 'include("Research/Research.jl"); using .Research; wp(; 
k=0.0, 
l=2*pi/30.0e3, 
m=-2*pi/1.0e3, 
branch=1, 
rx=0.0, 
ry=0.0, 
rz=15.0e3,
a0=0.7,
version=3.0,
coriolis_frequency=2.0e-3,
lx=900.0e3,
ly=30.0e3,
lz=80.0e3,
x0=0.0,
y0=0.0,
z0=30.0e3,
x_size=128,
y_size=32,
z_size=2400,
npx=16,
npy=4,
output_file="inertia_gw_wr_version3.h5", 
output_interval=100.0, 
tmax=100.0, 
turbulence_scheme=:TKEScheme, 
tracer_setup=:TracerOn,
prepare_restart=false,
momentum_coupling = true,
entropy_coupling = true,
turbulence_filter_order = 3,
restart = false,
input_file = "inertia_gw_wr_with_coupling.h5")' &> inertia_gw_wr.log

exit 0