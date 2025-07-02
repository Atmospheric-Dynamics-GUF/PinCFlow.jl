"""
    Boundaries

Module for setting boundary conditions on different variable types and field dimensions.
Handles periodic boundaries, solid wall boundaries, and MPI halo exchanges.
"""
module Boundaries

using ..Types
using ..MPIOperations

"""
Abstract type for boundary variable categories.
"""
abstract type AbstractBoundaryVariables end

"""
Boundary variables for predictand fields (rho, rhop, u, v, w, pip).
"""
struct BoundaryPredictands <: AbstractBoundaryVariables end

"""
Boundary variables for reconstruction fields.
"""
struct BoundaryReconstructions <: AbstractBoundaryVariables end

"""
Boundary variables for flux fields.
"""
struct BoundaryFluxes <: AbstractBoundaryVariables end

"""
Boundary variables for gravity wave integral fields.
"""
struct BoundaryGWIntegrals <: AbstractBoundaryVariables end

"""
Boundary variables for gravity wave tendency fields.
"""
struct BoundaryGWTendencies <: AbstractBoundaryVariables end

include("set_boundaries!.jl")
include("set_compressible_meridional_boundaries!.jl")
include("set_compressible_vertical_boundaries!.jl")
include("set_compressible_zonal_boundaries!.jl")
include("set_tracer_meridional_boundaries!.jl")
include("set_tracer_vertical_boundaries!.jl")
include("set_tracer_zonal_boundaries!.jl")
include("set_ice_meridional_boundaries!.jl")
include("set_ice_vertical_boundaries!.jl")
include("set_ice_zonal_boundaries!.jl")
include("set_turbulence_meridional_boundaries!.jl")
include("set_turbulence_vertical_boundaries!.jl")
include("set_turbulence_zonal_boundaries!.jl")
include("set_meridional_boundaries_of_field!.jl")
include("set_meridional_boundaries!.jl")
include("set_vertical_boundaries_of_field!.jl")
include("set_vertical_boundaries!.jl")
include("set_zonal_boundaries_of_field!.jl")
include("set_zonal_boundaries!.jl")

export BoundaryPredictands,
    BoundaryReconstructions,
    BoundaryFluxes,
    BoundaryGWIntegrals,
    BoundaryGWTendencies

export set_boundaries!,
    set_meridional_boundaries_of_field!,
    set_meridional_boundaries!,
    set_vertical_boundaries_of_field!,
    set_vertical_boundaries!,
    set_zonal_boundaries_of_field!,
    set_zonal_boundaries!

end
