module Boundaries

using ..Types
using ..MPIOperations

abstract type AbstractBoundaryVariables end

struct BoundaryPredictands <: AbstractBoundaryVariables end
struct BoundaryReconstructions <: AbstractBoundaryVariables end
struct BoundaryFluxes <: AbstractBoundaryVariables end
struct BoundaryGWIntegrals <: AbstractBoundaryVariables end
struct BoundaryGWTendencies <: AbstractBoundaryVariables end

include("set_boundaries!.jl")
include("set_compressible_meridional_boundaries!.jl")
include("set_compressible_vertical_boundaries!.jl")
include("set_compressible_zonal_boundaries!.jl")
include("set_meridional_boundaries_of_field!.jl")
include("set_meridional_boundaries_of_reduced_field!.jl")
include("set_meridional_boundaries!.jl")
include("set_vertical_boundaries_of_field!.jl")
include("set_vertical_boundaries_of_reduced_field!.jl")
include("set_vertical_boundaries!.jl")
include("set_zonal_boundaries_of_field!.jl")
include("set_zonal_boundaries_of_reduced_field!.jl")
include("set_zonal_boundaries!.jl")

export BoundaryPredictands,
    BoundaryReconstructions,
    BoundaryFluxes,
    BoundaryGWIntegrals,
    BoundaryGWTendencies

export set_boundaries!,
    set_meridional_boundaries_of_field!,
    set_meridional_boundaries_of_reduced_field!,
    set_meridional_boundaries!,
    set_vertical_boundaries_of_field!,
    set_vertical_boundaries_of_reduced_field!,
    set_vertical_boundaries!,
    set_zonal_boundaries_of_field!,
    set_zonal_boundaries_of_reduced_field!,
    set_zonal_boundaries!

end
