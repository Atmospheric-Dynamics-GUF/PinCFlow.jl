"""
```julia
Boundaries
```

Module for enforcing boundary conditions for different variable types and field dimensions.

Handles periodic boundaries in the horizontal and solid-wall boundaries in the vertical, as well as MPI communication in all dimensions.

# See also

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.MPIOperations`](@ref)
"""
module Boundaries

using ..Types
using ..MPIOperations
using ..PinCFlow

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
BoundaryWKBTendencies <: AbstractBoundaryVariables
```

Boundary-variable category for gravity-wave-tendency fields.
"""
struct BoundaryWKBTendencies <: AbstractBoundaryWKBVariables end

include("set_boundaries!.jl")
include("set_tracer_meridional_boundaries!.jl")
include("set_tracer_vertical_boundaries!.jl")
include("set_tracer_zonal_boundaries!.jl")
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
