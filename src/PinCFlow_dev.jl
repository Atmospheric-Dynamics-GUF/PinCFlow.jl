module PinCFlow_dev

include("types.jl")
include("matrix_solvers.jl")

include("semi_discretization.jl")

include("init.jl")
include("atmosphere.jl")
include("boundary.jl")
include("fluxes.jl")
include("update.jl")
include("time_step.jl")
export initialize_values, initialize_atmosphere!, initialize_variables!, setBoundary!,
    reconstruction!, compute_fluxes!

export time_discretization, massUpdate_rho!, massUpdate_rhop!    

export SemiDiscretization, pincflow

# debugging
export momentumPredictor_u!, momentumPredictor_v!, momentumPredictor_w!

end # module PinCFlow_dev
