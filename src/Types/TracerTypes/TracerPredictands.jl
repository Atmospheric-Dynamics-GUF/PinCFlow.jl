"""
```julia
TracerPredictands{A <: AbstractArray{<:AbstractFloat, 3}}
```

Arrays for tracers.

```julia
TracerPredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)::TracerPredictands
```

Construct a `TracerPredictands` instance with dimensions and initial values depending on the general configuration of tracer transport, by dispatching to the appropriate method.

```julia
TracerPredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    tracer_setup::NoTracer,
    variables::Variables,
)::TracerPredictands
```

Construct a `TracerPredictands` instance with zero-size arrays for configurations without tracer transport.

```julia
TracerPredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    tracer_setup::TracerOn,
    variables::Variables,
)::TracerPredictands
```

Construct a `TracerPredictands` instance with a tracer initialized by the function `initial_tracer` in `namelists.tracer`. The tracer field is multiplied by the density.

# Fields

  - `chi::A`: Non-dimensional tracer.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `atmosphere`: Atmospheric-background fields.

  - `grid`: Collection of parameters and fields describing the grid.

  - `tracer_setup`: General tracer-transport configuration.

  - `variables`: Container for arrays needed for the prediction of the prognostic variables.

# See also

  - [`PinCFlow.Types.FoundationalTypes.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.set_meridional_boundaries_of_field!`](@ref)
"""
struct TracerPredictands{A <: AbstractArray{<:AbstractFloat, 3}}
    chi::A
end

function TracerPredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)::TracerPredictands
    (; tracer_setup) = namelists.tracer

    return TracerPredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        tracer_setup,
        variables,
    )
end

function TracerPredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    tracer_setup::NoTracer,
    variables::Variables,
)::TracerPredictands
    return TracerPredictands(
        [zeros(0, 0, 0) for field in fieldnames(TracerPredictands)]...,
    )
end

function TracerPredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    tracer_setup::TracerOn,
    variables::Variables,
)::TracerPredictands
    (; nxx, nyy, nzz, i0, i1, j0, j1) = domain
    (; x, y, zc) = grid
    (; lref) = constants
    (; rhobar) = atmosphere
    (; rho) = variables.predictands
    (; lref) = constants
    (; initial_tracer) = namelists.tracer

    chi = zeros(nxx, nyy, nzz)
    @ivy for k in 1:nzz, j in j0:j1, i in i0:i1
        chi[i, j, k] =
            initial_tracer(x[i] * lref, y[j] * lref, zc[i, j, k] * lref)
    end
    set_zonal_boundaries_of_field!(chi, namelists, domain)
    set_meridional_boundaries_of_field!(chi, namelists, domain)

    chi .*= rho .+ rhobar

    return TracerPredictands(chi)
end
