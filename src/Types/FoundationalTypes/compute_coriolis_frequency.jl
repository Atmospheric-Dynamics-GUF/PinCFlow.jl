"""
```julia
compute_coriolis_frequency(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    coriolis_mode::FPlane,
)
```

Set the Coriolis parameter to ``f = f_0``, with ``f_0`` being given by `namelists.atmosphere.coriolis_frequency`.

# Arguments

  - `namelists`: Simulation parameters
  - `constants`: Simulation constants
  - `domain`: Computational domain
  - `grid`: not used
  - `coriolis_mode`: Type dispatch on mode for coriolis force computation

# Returns

  - `::AbstractVector{<:AbstractFloat}`: Coriolis parameter with no meridional dependence.
"""
function compute_coriolis_frequency(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    coriolis_mode::FPlane,
)
    (; coriolis_frequency) = namelists.atmosphere
    (; tref) = constants
    (; nyy) = domain
    return coriolis_frequency .* tref .* ones(nyy)
end
