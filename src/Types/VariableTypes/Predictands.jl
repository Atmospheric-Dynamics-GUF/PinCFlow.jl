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
)::Predictands
```

Construct a `Predictands` instance with dimensions and initial values depending on which test case is initialized, by dispatching to the appropriate method.

```julia
Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    testcase::AbstractTestCase,
)::Predictands
```

Construct a `Predictands` instance. The mass-weighted potential temperature `p` is constructed depending on the dynamic equations (see `set_p`).

The wind is initialized with ``\\boldsymbol{u}_0`` (given by `namelists.atmosphere.backgroundflow_dim`) everywhere, whereas the density fluctuations and Exner-pressure fluctuations are initialized with zero. The array for the mass-weighted potential temperature is constructed with size `(0, 0, 0)`.

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

  - `model`: Dynamic equations.

  - `testcase`: Test case on which the current simulation is based.
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
    (; testcase) = namelists.setting
    return Predictands(namelists, constants, domain, atmosphere, grid, testcase)
end

function Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    testcase::AbstractTestCase,
)::Predictands
    (; backgroundflow_dim) = namelists.atmosphere
    (; model) = namelists.setting
    (; uref) = constants
    (; nxx, nyy, nzz) = domain
    (; pstrattfc) = atmosphere

    # Initialize the predictands.
    (rho, rhop, u, v, w, pip) = (zeros(nxx, nyy, nzz) for i in 1:6)
    p = set_p(model, nxx, nyy, nzz, pstrattfc)

    # Set the initial winds.
    @ivy u .= backgroundflow_dim[1] ./ uref
    @ivy v .= backgroundflow_dim[2] ./ uref
    @ivy w .= backgroundflow_dim[3] ./ uref

    return Predictands(rho, rhop, u, v, w, pip, p)
end
