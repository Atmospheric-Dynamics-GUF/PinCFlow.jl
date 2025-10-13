"""
```julia
Predictands{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
}
```

Arrays for prognostic variables.

```julia
Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
)::Predictands
```

Construct a `Predictands` instance. The mass-weighted potential temperature `p` is constructed depending on the dynamic equations (see `set_p`).

The wind is initialized with ``\\boldsymbol{u}_0`` (given by `namelists.atmosphere.initial_wind`) everywhere, whereas the density fluctuations and Exner-pressure fluctuations are initialized with zero. The array for the mass-weighted potential temperature is constructed with size `(0, 0, 0)`.

# Fields

  - `rho::A`: Density.

  - `rhop::A`: Density-fluctuations.

  - `u::A`: Zonal wind.

  - `v::A`: Meridional wind.

  - `w::A`: Transformed vertical wind.

  - `pip::A`: Exner-pressure fluctuations.

  - `p::B`: Mass-weighted potential temperature.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `atmosphere`: Atmospheric-background fields.

  - `grid`: Collection of parameters and fields describing the grid.
"""
struct Predictands{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
}
    rho::A
    rhop::A
    u::A
    v::A
    w::A
    pip::A
    p::B
end

function Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
)::Predictands
    (; initial_wind) = namelists.atmosphere
    (; model) = namelists.setting
    (; uref) = constants
    (; nxx, nyy, nzz) = domain
    (; pbar) = atmosphere

    # Initialize the predictands.
    (rho, rhop, u, v, w, pip) = (zeros(nxx, nyy, nzz) for i in 1:6)
    p = set_p(model, nxx, nyy, nzz, pbar)

    # Set the initial winds.
    @ivy u .= initial_wind[1] ./ uref
    @ivy v .= initial_wind[2] ./ uref
    @ivy w .= initial_wind[3] ./ uref

    return Predictands(rho, rhop, u, v, w, pip, p)
end
