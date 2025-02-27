module PinCFlow_dev

include("types.jl")
include("matrix_solvers.jl")

include("semi_discretization.jl")

include("init.jl")
include("atmosphere.jl")
include("boundary.jl")
include("fluxes.jl")

export initialize_values, initialize_atmosphere!, initialize_variables!, setBoundary!,
    reconstruction!, compute_fluxes!

export SemiDiscretization

end # module PinCFlow_dev
