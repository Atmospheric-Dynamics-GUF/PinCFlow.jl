module PinCFlow

using OffsetArrays
using LinearAlgebra
using MPI

# using SimpleUnPack
# using TimerOutputs
# using TrixiBas

# Include types.
include("types/namelists/domain_namelist.jl")
include("types/namelists/output_namelist.jl")
include("types/namelists/setting_namelist.jl")
include("types/namelists/discretization_namelist.jl")
include("types/namelists/poisson_namelist.jl")
include("types/namelists/atmosphere_namelist.jl")
include("types/namelists/grid_namelist.jl")
include("types/namelists/sponge_namelist.jl")
include("types/namelists/boundaries_namelist.jl")
include("types/namelists/namelists.jl")
include("types/time.jl")
include("types/constants.jl")
include("types/domain.jl")
include("types/grid.jl")
include("types/atmosphere.jl")
include("types/sponge.jl")
include("types/poisson/tensor.jl")
include("types/poisson/operator.jl")
include("types/poisson/preconditioner.jl")
include("types/poisson/bicgstab.jl")
include("types/poisson/correction.jl")
include("types/poisson/poisson.jl")
include("types/variables/predictands.jl")
include("types/variables/tendencies.jl")
include("types/variables/backups.jl")
include("types/variables/auxiliaries.jl")
include("types/variables/reconstructions.jl")
include("types/variables/fluxes.jl")
include("types/variables/variables.jl")
include("types/state.jl")

# Include boundary functions.
include("boundaries/set_boundaries.jl")
include("boundaries/set_meridional_boundaries_of_field.jl")
include("boundaries/set_meridional_boundaries.jl")
include("boundaries/set_vertical_boundaries.jl")
include("boundaries/set_zonal_boundaries_of_field.jl")
include("boundaries/set_zonal_boundaries.jl")

# Include flux functionss.
include("fluxes/apply_1d_muscl.jl")
include("fluxes/apply_3d_muscl.jl")
include("fluxes/compute_flux.jl")
include("fluxes/compute_fluxes.jl")
include("fluxes/reconstruct.jl")

# Include MPI functions.
include("mpi/compute_global_dot_product.jl")
include("mpi/set_meridional_halos_of_field.jl")
include("mpi/set_zonal_halos_of_field.jl")

# Include update functions.
include("update/apply_unified_sponge.jl")
include("update/compute_sponge.jl")
include("update/compute_stress_tensor.jl")
include("update/compute_time_step.jl")
include("update/compute_vertical_wind.jl")
include("update/transform.jl")
include("update/update.jl")

# Include Poisson functions.
include("poisson/apply_operator.jl")
include("poisson/apply_preconditioner.jl")
include("poisson/apply_bicgstab.jl")
include("poisson/compute_operator.jl")
include("poisson/compute_rhs.jl")
include("poisson/solve_poisson.jl")
include("poisson/correct.jl")
include("poisson/apply_corrector.jl")

# Include integration functions.
include("integration/integrate.jl")

# pincflow_test_dir() = joinpath(dirname(pathof(PinCFlow)), "..", "test")
# pincflow_examples_dir() = joinpath(dirname(pathof(PinCFlow)), "..", "examples")

# export pincflow, Corrector
# export pincflow_test_dir, pincflow_examples_dir
# export Grid, Constants, Atmosphere, Domain, Variables, Fluxes, State
# export Corrector
#
export Namelists, State

end
