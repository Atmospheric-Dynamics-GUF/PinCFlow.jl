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
using ..PinCFlow

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

"""
```julia
X
```

Singleton for ``\\widehat{x}``-axis along which a calculation should be performed.
"""
struct X end

"""
```julia
Y
```

Singleton for ``\\widehat{y}``-axis along which a calculation should be performed.
"""
struct Y end

"""
```julia
Z
```

Singleton for ``\\widehat{z}``-axis along which a calculation should be performed.
"""
struct Z end

include("apply_lhs_sponge!.jl")
include("compute_buoyancy_factor.jl")
include("compute_compressible_wind_factor.jl")
include("compute_pressure_gradient.jl")
include("compute_sponges!.jl")
include("compute_stress_tensor.jl")
include("compute_vertical_wind.jl")
include("compute_volume_force.jl")
include("transform.jl")
include("update!.jl")
include("conductive_heating.jl")
include("compute_momentum_diffusion_terms.jl")

export LHS, RHS, X, Y, Z

export apply_lhs_sponge!,
    compute_buoyancy_factor,
    compute_compressible_wind_factor,
    compute_pressure_gradient,
    compute_sponges!,
    compute_stress_tensor,
    compute_vertical_wind,
    compute_volume_force,
    conductive_heating,
    compute_momentum_diffusion_terms,
    update!

end
