"""
```julia
Update
```

Module for integrating the prognostic equations.

Provides functions for updating the prognostic variables at the various stages of the semi-implicit time scheme.

# See also

- [`PinCFlow.Types`](@ref)
- [`PinCFlow.Boundaries`](@ref)
"""
module Update

using MPI
using ..Types
using ..Boundaries

"""
```julia
Cartesian
```

Singleton for transformations from the terrain-following system to the Cartesian one.
"""
struct Cartesian end

"""
```julia
Transformed
```

Singleton for transformations from the Cartesian system to the terrain-following one.
"""
struct Transformed end

"""
```julia
LHS
```

Singleton for the integration of the left-hand side of an equation.
"""
struct LHS end

"""
```julia
RHS
```

Singleton for the integration of the right-hand side of an equation.
"""
struct RHS end

include("apply_unified_sponge!.jl")
include("compute_compressible_buoyancy_factor.jl")
include("compute_compressible_wind_factor.jl")
include("compute_sponge!.jl")
include("compute_stress_tensor.jl")
include("compute_vertical_wind.jl")
include("compute_volume_force.jl")
include("transform.jl")
include("update!.jl")

export LHS, RHS

export apply_unified_sponge!,
    compute_compressible_buoyancy_factor,
    compute_compressible_wind_factor,
    compute_sponge!,
    compute_stress_tensor,
    compute_vertical_wind,
    compute_volume_force,
    update!

end
