module PinCFlow

using OffsetArrays
using LinearAlgebra
using NetCDF
using MPI

using SimpleUnPack
using TimerOutputs
using TrixiBase

# Namelists
include("namelists/atmosphere_namelist.jl")
include("namelists/boundaries_namelist.jl")
include("namelists/discretization_namelist.jl")
include("namelists/domain_namelist.jl")
include("namelists/grid_namelist.jl")
include("namelists/namelists.jl")
include("namelists/output_namelist.jl")
include("namelists/poisson_namelist.jl")
include("namelists/settings_namelist.jl")
include("namelists/sponge_namelist.jl")

# Atmosphere
include("atmosphere/atmosphere.jl")
include("atmosphere/constants.jl")
include("atmosphere/grid.jl")
include("atmosphere/sponge.jl")

# Boundaries
include("boundaries/boundaries.jl")
include("boundaries/meridional_boundaries.jl")
include("boundaries/vertical_boundaries.jl")
include("boundaries/zonal_boundaries.jl")

# Fluxes
include("fluxes/mass_fluxes.jl")
include("fluxes/momentum_fluxes.jl")
include("fluxes/reconstruction.jl")

# Integration
include("integration/integration.jl")
include("integration/state.jl")

# Output
include("output/output.jl")

# Poisson
include("poisson/poisson.jl")

# Update
include("update/sponge.jl")
include("update/time.jl")
include("update/time_loop.jl")
include("update/time_step.jl")
include("update/update.jl")

# Variables
include("variables/stress_tensor.jl")
include("variables/variables.jl")
include("variables/vertical_wind.jl")

pincflow_test_dir() = joinpath(dirname(pathof(PinCFlow)), "..", "test")
pincflow_examples_dir() = joinpath(dirname(pathof(PinCFlow)), "..", "examples")

export pincflow, Corrector
export pincflow_test_dir, pincflow_examples_dir
export Grid, Constants, Atmosphere, Domain, Variables, Fluxes, State
export Corrector

end
