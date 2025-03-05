module PinCFlow_dev

include("types.jl")
include("matrix_solvers.jl")

include("semi_discretization.jl")

include("init.jl")
include("sponge.jl")
include("sponge_new.jl")
include("atmosphere.jl")
include("boundary.jl")
include("fluxes.jl")
include("time_step.jl")
include("poisson.jl")
include("update.jl")

pincflow_test_dir() = joinpath(dirname(pathof(PinCFlow_dev)), "..", "test")

export initialize_values,
    initialize_atmosphere!,
    initialize_variables!,
    setBoundary!,
    reconstruction!,
    compute_fluxes!,
    vertWind,
    time_loop!

export time_discretization, massUpdate_rho!, massUpdate_rhop!

export SemiDiscretization, pincflow, Corrector

export pincflow_test_dir

# debugging
export Corrector

end # module PinCFlow_dev
