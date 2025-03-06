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

include("namelist_interface.jl")

pincflow_test_dir() = joinpath(dirname(pathof(PinCFlow_dev)), "..", "test")
pincflow_examples_dir() = joinpath(dirname(pathof(PinCFlow_dev)), "..", "examples")

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

export setup_semidiscretization, setup_atmosphere_list, setup_grid_list, setup_output_list,
       setup_debugging_list, setup_test_case_list, setup_model, setup_solver_list,
       setup_poisson_solver_list, setup_topography_list, setup_boundary_list

export pincflow_test_dir, pincflow_examples_dir

# debugging
export Corrector

end # module PinCFlow_dev
