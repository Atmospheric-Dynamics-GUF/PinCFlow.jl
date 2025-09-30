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
    tracersetup::AbstractTracer,
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
    tracersetup::LinearTracer,
    variables::Variables,
)::TracerPredictands
```

Construct a `TracerPredictands` instance with an initialized non-dimensional tracer linearly increasing with altitude. The tracer field is multiplied by the density.

# Fields

  - `chi::A`: Non-dimensional tracer.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `atmosphere`: Atmospheric-background fields.

  - `grid`: Collection of parameters and fields describing the grid.

  - `tracersetup`: General tracer-transport configuration.

  - `variables`: Container for arrays needed for the prediction of the prognostic variables.
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
    (; tracersetup) = namelists.tracer

    return TracerPredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        tracersetup,
        variables,
    )
end

function TracerPredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    tracersetup::AbstractTracer,
    variables::Variables,
)::TracerPredictands
    chi = zeros(0, 0, 0)

    return TracerPredictands(chi)
end

function TracerPredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    tracersetup::LinearTracer,
    variables::Variables,
)::TracerPredictands
    (; nxx, nyy, nzz) = domain
    (; zc) = grid
    (; lref) = constants
    (; rhobar) = atmosphere
    (; rho) = variables.predictands
    (; lref) = constants
    (; alphatracer) = namelists.tracer

    chi = zeros(nxx, nyy, nzz)
    chi .= alphatracer .* lref .* zc

    chi .= chi .* (rho .+ rhobar)

    return TracerPredictands(chi)
end
