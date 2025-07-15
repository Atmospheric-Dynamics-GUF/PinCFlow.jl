"""
```julia
compute_coriolis_parameter(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    coriolis_mode::FPlane,
)
```

Set the Coriolis parameter to a constant.

# Arguments

  - `namelists`: Simulation parameters
  - `constants`: Simulation constants
  - `domain`: Computational domain
  - `grid`: not used
  - `coriolis_mode`: Type dispatch on mode for coriolis force computation

# Returns

  - Array with size `nyy`, initialized with coriolis_frequency * t_ref
"""
function compute_coriolis_parameter(
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
