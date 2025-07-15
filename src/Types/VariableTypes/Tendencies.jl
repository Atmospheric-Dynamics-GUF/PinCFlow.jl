"""
```julia
Tendencies{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
}
```

Time tendency storage for the prognostic variables.

Holds the time derivatives (tendencies) computed during each Runge-Kutta stage
for all prognostic fields. These tendencies are accumulated from various
processes including advection, and pressure forces.

# Fields

  - `drho::A`: Total density tendency [∂ρ/∂t]
  - `drhop::A`: Density perturbation tendency [∂ρ'/∂t]
  - `du::A`: Zonal velocity tendency [∂u/∂t]
  - `dv::A`: Meridional velocity tendency [∂v/∂t]
  - `dw::A`: Vertical velocity tendency [∂w/∂t]
  - `dpip::A`: Pressure correction tendency [∂π'/∂t]
  - `dp::B`: Pressure tendency (compressible model only) [∂p/∂t]

# Constructors

    Tendencies(namelists::Namelists, domain::Domain)
    Tendencies(domain::Domain, model::AbstractModel)
    Tendencies(domain::Domain, model::Compressible)

# Model Dependencies

## AbstractModel (Boussinesq, PseudoIncompressible)

  - Uses 6 tendency fields (dp unused, set to empty array)
  - Pressure handled through dpip correction field
  - Density split into background + perturbation components

## Compressible

  - Uses all 7 tendency fields including dp
  - Explicit pressure evolution equation
"""
struct Tendencies{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
}
    drho::A
    drhop::A
    du::A
    dv::A
    dw::A
    dpip::A
    dp::B
end

"""
```julia
Tendencies(namelists::Namelists, domain::Domain)
```

Create tendencies storage from configuration.

Dispatches to model-specific constructor based on equation set specified
in namelists.setting.model.

# Arguments

  - `namelists`: Configuration including model type
  - `domain`: Domain decomposition for array sizing

# Returns

  - `Tendencies`: Initialized tendency storage for the specified model
"""
function Tendencies(namelists::Namelists, domain::Domain)
    (; model) = namelists.setting
    return Tendencies(domain, model)
end

"""
```julia
Tendencies(domain::Domain, model::AbstractModel)
```

Create tendencies for non-compressible models (Boussinesq, PseudoIncompressible).

# Arguments

  - `domain`: Domain information for array dimensions
  - `model`: Model type (excludes Compressible)

# Returns

  - `Tendencies`: Storage with 6 active fields, dp field unused (empty array)
"""
function Tendencies(domain::Domain, model::AbstractModel)
    (; nxx, nyy, nzz) = domain

    # Initialize the tendencies.
    (drho, drhop, du, dv, dw, dpip) = (zeros(nxx, nyy, nzz) for i in 1:6)
    dp = zeros(0, 0, 0)

    # Return a Variables instance.
    return Tendencies(drho, drhop, du, dv, dw, dpip, dp)
end

"""
```julia
Tendencies(domain::Domain, model::Compressible)
```

Create tendencies for fully compressible model.

# Arguments

  - `domain`: Domain information for array dimensions
  - `model`: Compressible model type

# Returns

  - `Tendencies`: Storage with all 7 fields active including explicit pressure tendency
"""
function Tendencies(domain::Domain, model::Compressible)
    (; nxx, nyy, nzz) = domain

    # Initialize the tendencies.
    (drho, drhop, du, dv, dw, dpip, dp) = (zeros(nxx, nyy, nzz) for i in 1:7)

    # Return a Variables instance.
    return Tendencies(drho, drhop, du, dv, dw, dpip, dp)
end
