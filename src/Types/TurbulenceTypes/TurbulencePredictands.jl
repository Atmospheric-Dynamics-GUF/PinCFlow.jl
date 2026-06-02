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

Construct a `TurbulencePredictands` instance with dimensions and initial values depending on the general turbulence parameterization configuration, by dispatching to the appropriate method.

```julia
TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulence_scheme::Val{:NoTurbulence},
    variables::Variables,
)::TurbulencePredictands
```

Construct a `TurbulencePredictands` instance with zero-size arrays for configurations without turbulence parameterization.

```julia 
TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulence_scheme::Val{:TKEScheme},
    variables::Variables,
)::TurbulencePredictands
```

Construct a `TurbulencePredictands` instance with the turbulent kinetic energy initialized by the function `initial_tke` in `namelists.turbulence`. The turbulence field is multiplied by the density.

# Fields

  - `tke::A`: Non-dimensional turbulent kinetic energy.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `atmosphere`: Atmospheric-background fields.

  - `grid`: Collection of parameters and fields describing the grid.

  - `turbulence_scheme`: General turbulence parameterization configuration.

  - `variables`: Container for arrays needed for the prediction of the prognostic variables.

# See also

  - [`PinCFlow.Types.FoundationalTypes.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.set_meridional_boundaries_of_field!`](@ref)
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

    @dispatch_turbulence_scheme return TurbulencePredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        Val(turbulence_scheme),
        variables,
    )
end

function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulence_scheme::Val{:NoTurbulence},
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
    turbulence_scheme::Val{:TKEScheme},
    variables::Variables,
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
