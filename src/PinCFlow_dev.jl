module PinCFlow_dev
using TrixiBase # get the `timer` object

# Namelists
include("namelists/atmosphere_namelist.jl")
include("namelists/boundary_namelist.jl")
include("namelists/discretization_namelist.jl")
include("namelists/domain_namelist.jl")
include("namelists/grid_namelist.jl")
include("namelists/namelists.jl")
include("namelists/output_namelist.jl")
include("namelists/poisson_namelist.jl")
include("namelists/settings_namelist.jl")
include("namelists/sponge_namelist.jl")
#
# Atmosphere
include("atmosphere/atmosphere.jl")
include("atmosphere/constants.jl")
include("atmosphere/grid.jl")
include("atmosphere/sponge.jl")

# Boundaries
include("boundaries/boundary.jl")
include("boundaries/domain.jl")

# Fluxes
include("fluxes/calc_flux.jl")
include("fluxes/reconstruction.jl")

# IO
include("io/save_solution.jl")

# Poisson
include("poisson/poisson.jl")

# Update
include("update/time.jl")
include("update/time_loop.jl")
include("update/time_step.jl")
include("update/update.jl")

# Variables
include("variables/variables.jl")

# Model
include("model.jl")

pincflow_test_dir() = joinpath(dirname(pathof(PinCFlow_dev)), "..", "test")
pincflow_examples_dir() = joinpath(dirname(pathof(PinCFlow_dev)), "..", "examples")


export pincflow, Corrector


export pincflow_test_dir, pincflow_examples_dir


export Grid, Constants, Atmosphere, Domain, Variables, Fluxes, State

# debugging
export Corrector

end # module PinCFlow_dev
