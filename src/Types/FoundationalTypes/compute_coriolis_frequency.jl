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

- `namelists`: Namelists with all model parameters.
- `constants`: Physical constants and reference values.
- `domain`: Collection of domain-decomposition and MPI-communication parameters.
- `grid`: Collection of parameters and fields describing the grid.
- `coriolis_mode`: Approximation used for the Coriolis frequency.

# Returns

- `::AbstractVector{<:AbstractFloat}`: Coriolis parameter with the meridional dependence specified by `coriolis_mode`.
"""
function compute_coriolis_frequency end

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
