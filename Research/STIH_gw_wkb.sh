#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=STIHwkb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --hint=nomultithread
#SBATCH --mem=90G
#SBATCH --time=0-00:10:00
#SBATCH --mail-type=ALL
#SBATCH --account=bb1097

set -x

# Set Intel MPI configuration on compute partition.
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
srun --distribution=block:cyclic julia --project=Research -e 'include("Research/Research.jl"); using .Research; wkb_wp(; 
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
x_size=1,
y_size=1,
z_size=100,
output_file="STIH_gw_wkb_lowres.h5", 
output_interval = 360.0, 
tmax = 3600.0 * 15, 
turbulence_scheme = :NoTurbulence, 
model = :PseudoIncompressible,
tracer_setup = :NoTracer,
next_order_impact = true,
turbulence_impact = true,
leading_order_impact = true,
turbulence_filter_order = 0,
momentum_coupling = true,
entropy_coupling = true)' &> STIH_gw_wkb.log

exit 0