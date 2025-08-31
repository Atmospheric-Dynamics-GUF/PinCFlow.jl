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

Construct a `Predictands` instance with dimensions and initial values depending on whether or not the model is compressible and which test case is initialized, by dispatching to the appropriate method.

```julia
Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    model::AbstractModel,
    testcase::AbstractTestCase,
)::Predictands
```

Construct a `Predictands` instance in non-compressible modes and non-wave-packet test cases.

The wind is initialized with ``\\boldsymbol{u}_0`` (given by `namelists.atmosphere.backgroundflow_dim`) everywhere, whereas the density fluctuations and Exner-pressure fluctuations are initialized with zero. The array for the mass-weighted potential temperature is constructed with size `(0, 0, 0)`.

```julia
Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    model::Compressible,
    testcase::AbstractTestCase,
)::Predictands
```

Construct a `Predictands` instance in compressible mode and non-wave-packet test cases.

The wind is initialized with ``\\boldsymbol{u}_0`` (given by `namelists.atmosphere.backgroundflow_dim`) everywhere, whereas the density fluctuations and Exner-pressure fluctuations are initialized with zero. The mass-weighted potential temperature is initialized with the corresponding array in `atmosphere`.

```julia
Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    model::PseudoIncompressible,
)::Predictands
```

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
    (; model, testcase) = namelists.setting
    return Predictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        model,
        testcase,
    )
end

function Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    model::AbstractModel,
    testcase::AbstractTestCase,
)::Predictands
    (; backgroundflow_dim) = namelists.atmosphere
    (; uref) = constants
    (; nxx, nyy, nzz) = domain

    # Initialize the predictands.
    (rho, rhop, u, v, w, pip) = (zeros(nxx, nyy, nzz) for i in 1:6)
    p = zeros(0, 0, 0)

    # Set the initial winds.
    u .= backgroundflow_dim[1] / uref
    v .= backgroundflow_dim[2] / uref
    w .= backgroundflow_dim[3] / uref

    # Return a Predictands instance.
    return Predictands(rho, rhop, u, v, w, pip, p)
end

function Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    model::Compressible,
    testcase::AbstractTestCase,
)::Predictands
    (; backgroundflow_dim) = namelists.atmosphere
    (; uref) = constants
    (; nxx, nyy, nzz) = domain
    (; pstrattfc) = atmosphere

    # Initialize the predictands.
    (rho, rhop, u, v, w, pip, p) = (zeros(nxx, nyy, nzz) for i in 1:7)

    # Set the initial winds.
    u .= backgroundflow_dim[1] / uref
    v .= backgroundflow_dim[2] / uref
    w .= backgroundflow_dim[3] / uref

    # Set the initial mass-weighted potential temperature.
    p .= pstrattfc

    # Return a Predictands instance.
    return Predictands(rho, rhop, u, v, w, pip, p)
end
