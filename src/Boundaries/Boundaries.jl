"""
```julia
Boundaries
```

Module for enforcing boundary conditions for different variable types and field dimensions.

Handles periodic boundaries in the horizontal and solid-wall boundaries in the vertical, as well as MPI communication in all dimensions.

# Scopes

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.MPIOperations`](@ref)

  - [`PinCFlow.Types`](@ref)

# Exports

  - [`PinCFlow.Boundaries.BoundaryFluxes`](@ref)

  - [`PinCFlow.Boundaries.BoundaryPredictands`](@ref)

  - [`PinCFlow.Boundaries.BoundaryReconstructions`](@ref)

  - [`PinCFlow.Boundaries.BoundaryWKBIntegrals`](@ref)

  - [`PinCFlow.Boundaries.BoundaryWKBTendencies`](@ref)

  - [`PinCFlow.Boundaries.set_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_meridional_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)

  - [`PinCFlow.Boundaries.set_vertical_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)

  - [`PinCFlow.Boundaries.set_zonal_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
"""
module Boundaries

using ..Macros
using ..Types
using ..MPIOperations

"""
```julia
AbstractBoundaryVariables
```

Abstract type for boundary-variable categories.
"""
abstract type AbstractBoundaryVariables end

"""
```julia
AbstractBoundaryWKBVariables <: AbstractBoundaryVariables
```
"""
abstract type AbstractBoundaryWKBVariables <: AbstractBoundaryVariables end

"""
```julia
BoundaryPredictands <: AbstractBoundaryVariables
```

Boundary-variable category for predictand fields.
"""
struct BoundaryPredictands <: AbstractBoundaryVariables end

"""
```julia
BoundaryReconstructions <: AbstractBoundaryVariables
```

Boundary-variable category for reconstruction fields.
"""
struct BoundaryReconstructions <: AbstractBoundaryVariables end

"""
```julia
BoundaryFluxes <: AbstractBoundaryVariables
```

Boundary-variable category for flux fields.
"""
struct BoundaryFluxes <: AbstractBoundaryVariables end

"""
```julia
BoundaryWKBIntegrals <: AbstractBoundaryWKBVariables
```

Boundary-variable category for gravity-wave-integral fields.
"""
struct BoundaryWKBIntegrals <: AbstractBoundaryWKBVariables end

"""
```julia
BoundaryWKBTendencies <: AbstractBoundaryWKBVariables
```

Boundary-variable category for gravity-wave-tendency fields.
"""
struct BoundaryWKBTendencies <: AbstractBoundaryWKBVariables end

include("set_boundaries!.jl")
include("set_tracer_meridional_boundaries!.jl")
include("set_tracer_vertical_boundaries!.jl")
include("set_tracer_zonal_boundaries!.jl")
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
    BoundaryWKBIntegrals,
    BoundaryWKBTendencies

export set_boundaries!,
    set_meridional_boundaries_of_field!,
    set_meridional_boundaries!,
    set_vertical_boundaries_of_field!,
    set_vertical_boundaries!,
    set_zonal_boundaries_of_field!,
    set_zonal_boundaries!

end
