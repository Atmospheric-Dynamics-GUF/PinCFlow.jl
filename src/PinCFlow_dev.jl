module PinCFlow_dev
using TrixiBase # get the `timer` object

include("parameters.jl")
include("types.jl")
include("matrix_solvers.jl")
include("domain.jl")
include("constants.jl")
include("grid.jl")
include("variables.jl")
include("atmosphere.jl")
include("fluxes.jl")
include("poisson.jl")
include("model.jl")

include("boundary.jl")
include("semi_discretization.jl")

include("init.jl")
include("sponge.jl")
include("sponge_new.jl")
include("time_step.jl")
include("update.jl")

include("namelist_interface.jl")

pincflow_test_dir() = joinpath(dirname(pathof(PinCFlow_dev)), "..", "test")
pincflow_examples_dir() = joinpath(dirname(pathof(PinCFlow_dev)), "..", "examples")

export initialize_values,
    initialize_atmosphere!,
    initialize_variables!,
    initialize!,
    setBoundary!,
    reconstruction!,
    compute_fluxes!,
    vertWind,
    time_loop!

export time_discretization, massUpdate_rho!, massUpdate_rhop!

export SemiDiscretization, pincflow, Corrector

export setup_semidiscretization

export pincflow_test_dir, pincflow_examples_dir

export DomainParameters, OutputParameters, BoundaryParameters, TopographyParameters,
    AtmosphereParameters, PoissonSolverParameters, DiscretizationParameters,
    TestCaseParameters, DebugParameters, ModelParameters, Parameters

export Grid, Constants, Atmosphere, Domain, PoissonOperator, Variables, Fluxes, Model

# debugging
export Corrector

end # module PinCFlow_dev
