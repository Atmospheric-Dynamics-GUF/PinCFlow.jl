module MPIOperations

using MPI
using LinearAlgebra
using ..Types

include("compute_global_dot_product.jl")
include("set_meridional_halos_of_field!.jl")
include("set_meridional_halos_of_reduced_field!.jl")
include("set_zonal_halos_of_field!.jl")
include("set_zonal_halos_of_reduced_field!.jl")

export compute_global_dot_product,
    set_meridional_halos_of_field!,
    set_meridional_halos_of_reduced_field!,
    set_zonal_halos_of_field!,
    set_zonal_halos_of_reduced_field!

end