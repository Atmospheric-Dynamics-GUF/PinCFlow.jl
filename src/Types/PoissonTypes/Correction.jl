"""
```julia
Correction{A <: AbstractArray{<:AbstractFloat, 3}}
```

Storage for pressure correction terms applied to velocity fields.

# Fields

  - `corx::A`: Pressure correction in x-direction (nxx x nyy x nzz)
  - `cory::A`: Pressure correction in y-direction (nxx x nyy x nzz)

# Usage

Correction terms are computed from the pressure solution and applied by
[`PinCFlow.PoissonSolver.correct!`](@ref) to ensure mass conservation
in the velocity field.
"""
struct Correction{A <: AbstractArray{<:AbstractFloat, 3}}
    corx::A
    cory::A
end

"""
```julia
Correction(domain::Domain)
```

Initialize correction arrays sized according to extended domain.

# Arguments

  - `domain`: Domain specification with extended dimensions

# Returns

  - `::Correction`: `Correction` instance with zero-initialized arrays.
"""
function Correction(domain::Domain)

    # Get all necessary fields.
    (; nxx, nyy, nzz) = domain

    # Return a Correction instance.
    return Correction([zeros(nxx, nyy, nzz) for i in 1:2]...)
end
