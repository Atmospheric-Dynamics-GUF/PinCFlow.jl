"""
```julia
TurbulencePredictands{A <: AbstractArray{<:AbstractFloat, 3}}
```

Arrays for prognostic turbulence variables.

```julia
TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)::TurbulencePredictands
```

Construct a `TurbulencePredictands` instance with dimensions and initial values depending on the general configuration of turbulence physics, by dispatching to the appropriate method.

```julia
TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulence_scheme::NoTurbulence,
    variables::Variables,
)::TurbulencePredictands
```

Construct a `TurbulencePredictands` instance with zero-size arrays for configurations without turbulence physics.

```julia
TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulence_scheme::TKEScheme,
    variables::Variables,
)::TurbulencePredictands
```

Construct a `TurbulencePredictands` instance with both arrays initialized as ``t_\\mathrm{ref}^2 / \\left(10 L_\\mathrm{ref}^2\\right) \\rho`` (non-dimensionalized), where ``t_\\mathrm{ref}`` and ``L_\\mathrm{ref}`` are given by the properties `tref` and `lref` of `constants`, respectively.

# Fields

  - `tke::A`: Turbulent kinetic energy.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `atmosphere`: Atmospheric-background fields.

  - `grid`: Collection of parameters and fields describing the grid.

  - `turbulence_scheme`: General turbulence-physics configuration.

  - `variables`: Container for arrays needed for the prediction of the prognostic variables.
"""
struct TurbulencePredictands{A <: AbstractArray{<:AbstractFloat, 3}}
    tke::A
end

function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)::TurbulencePredictands
    (; turbulence_scheme) = namelists.turbulence

    return TurbulencePredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        turbulence_scheme,
        variables,
    )
end

function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulence_scheme::NoTurbulence,
    variables::Variables,
)::TurbulencePredictands
    tke = zeros(0, 0, 0)

    return TurbulencePredictands(tke)
end

function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulence_scheme::TKEScheme,
    variables::Variables,
)::TurbulencePredictands
    (; model) = namelists.atmosphere

    return TurbulencePredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        turbulence_scheme,
        variables,
        model,
    )
end

function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulence_scheme::TKEScheme,
    variables::Variables,
    model::Union{PseudoIncompressible, Compressible},
)::TurbulencePredictands
    (; i0, i1, j0, j1, nxx, nyy, nzz) = domain
    (; rhobar) = atmosphere
    (; x, y, zc) = grid
    (; rho) = variables.predictands
    (; lref, tref) = constants
    (; initial_tke) = namelists.turbulence

    tke = zeros(nxx, nyy, nzz)

    @ivy for k in 1:nzz, j in j0:j1, i in i0:i1
        xdim = x[i] * lref
        ydim = y[j] * lref
        zcdim = zc[i, j, k] * lref

        tke[i, j, k] =
            initial_tke(xdim, ydim, zcdim) * (tref^2.0) / (lref^2.0) *
            (rho[i, j, k] + rhobar[i, j, k])
    end

    for f! in
        (set_zonal_boundaries_of_field!, set_meridional_boundaries_of_field!)
        f!(tke, namelists, domain)
    end

    return TurbulencePredictands(tke)
end

function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulence_scheme::TKEScheme,
    variables::Variables,
    model::Boussinesq,
)::TurbulencePredictands
    (; i0, i1, j0, j1, nxx, nyy, nzz) = domain
    (; rhobar) = atmosphere
    (; x, y, zc) = grid
    (; rhop) = variables.predictands
    (; lref, tref) = constants
    (; initial_tke) = namelists.turbulence

    tke = zeros(nxx, nyy, nzz)

    @ivy for k in 1:nzz, j in j0:j1, i in i0:i1
        xdim = x[i] * lref
        ydim = y[j] * lref
        zcdim = zc[i, j, k] * lref

        tke[i, j, k] =
            initial_tke(xdim, ydim, zcdim) * (tref^2.0) / (lref^2.0) *
            (rhop[i, j, k] + rhobar[i, j, k])
    end

    for f! in
        (set_zonal_boundaries_of_field!, set_meridional_boundaries_of_field!)
        f!(tke, namelists, domain)
    end

    return TurbulencePredictands(tke)
end
