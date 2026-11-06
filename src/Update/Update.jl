"""
```julia
Update
```

Module for integrating the prognostic equations.

Provides functions for updating the prognostic variables at the various stages of the semi-implicit time scheme.

# Scopes

  - [`PinCFlow.Boundaries`](@ref)

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.Types`](@ref)

# Exports

  - [`PinCFlow.Update.KEK`](@ref)

  - [`PinCFlow.Update.KH`](@ref)

  - [`PinCFlow.Update.KM`](@ref)

  - [`PinCFlow.Update.LHS`](@ref)

  - [`PinCFlow.Update.RHS`](@ref)

  - [`PinCFlow.Update.apply_lhs_sponge!`](@ref)

  - [`PinCFlow.Update.check_tke!`](@ref)

  - [`PinCFlow.Update.compute_buoyancy_factor`](@ref)

  - [`PinCFlow.Update.compute_compressible_wind_factor`](@ref)

  - [`PinCFlow.Update.compute_momentum_diffusion_terms`](@ref)

  - [`PinCFlow.Update.compute_pressure_gradient`](@ref)

  - [`PinCFlow.Update.compute_sponges!`](@ref)

  - [`PinCFlow.Update.compute_stress_tensor!`](@ref)

  - [`PinCFlow.Update.compute_vertical_wind`](@ref)

  - [`PinCFlow.Update.compute_volume_force`](@ref)

  - [`PinCFlow.Update.conductive_heating`](@ref)

  - [`PinCFlow.Update.reset_thomas!`](@ref)

  - [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - [`PinCFlow.Update.turbulence_diffusion_coefficient`](@ref)

  - [`PinCFlow.Update.turbulent_diffusion!`](@ref)

  - [`PinCFlow.Update.update!`](@ref)
"""
module Update

using MPI
using ..Macros
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

"""
```julia
KM
```

Singleton for the eddy diffusion coefficient for momentum.
"""
struct KM end

"""
```julia
KH
```

Singleton for the eddy diffusion coefficient for heat.
"""
struct KH end

"""
```julia
KEK
```

Singleton for the eddy diffusion coefficient for turbulent kinetic energy.
"""
struct KEK end

include("check_tke!.jl")
include("turbulence_diffusion_coefficient.jl")
include("apply_lhs_sponge!.jl")
include("compute_buoyancy_factor.jl")
include("compute_compressible_wind_factor.jl")
include("compute_pressure_gradient.jl")
include("compute_sponges!.jl")
include("compute_stress_tensor!.jl")
include("compute_vertical_wind.jl")
include("compute_volume_force.jl")
include("transform.jl")
include("update!.jl")
include("conductive_heating.jl")
include("compute_momentum_diffusion_terms.jl")
include("thomas_algorithm!.jl")
include("reset_thomas!.jl")
include("turbulent_diffusion!.jl")

export LHS, RHS, KM, KH, KEK

export apply_lhs_sponge!,
    compute_buoyancy_factor,
    compute_compressible_wind_factor,
    compute_pressure_gradient,
    compute_sponges!,
    compute_stress_tensor!,
    compute_vertical_wind,
    compute_volume_force,
    conductive_heating,
    compute_momentum_diffusion_terms,
    thomas_algorithm!,
    reset_thomas!,
    turbulent_diffusion!,
    turbulence_diffusion_coefficient,
    check_tke!,
    update!

end
