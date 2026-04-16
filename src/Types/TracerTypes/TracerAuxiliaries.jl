"""
```julia
TracerAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```

Initial states of the tracers.

```julia
TracerAuxiliaries(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)::TracerAuxiliaries
```

Construct a `TracerAuxiliaries` instance with dimensions and initial values depending on the general configuration of tracer transport, by dispatching to the appropriate method.

```julia 
TracerAuxiliaries(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    tracer_setup::Val{:NoTracer},
    variables::Variables,
)::TracerAuxiliaries
```

Construct a `TracerAuxiliaries` instance with zero-size arrays for configurations without tracer transport.

```julia 
TracerAuxiliaries(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    tracer_setup::Val{:TracerOn},
    variables::Variables,
)::TracerAuxiliaries
```

Construct a `TracerAuxiliaries` instance with a background tracer initialized by the function `background_tracer` in `namelists.tracer`. The background tracer field is multiplied by the background density.

# Fields

  - `backgroundtracer::A`: Non-dimensional tracer to which the tracer distribution should be relaxed to if `namelists.tracer.apply_sponge_to_tracer` is set to `true`.

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
struct TracerAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    backgroundtracer::A
end

function TracerAuxiliaries(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)::TracerAuxiliaries
    (; tracer_setup) = namelists.tracer

    @dispatch_tracer_setup return TracerAuxiliaries(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        Val(tracer_setup),
        variables,
    )
end
function TracerAuxiliaries(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    tracer_setup::Val{:NoTracer},
    variables::Variables,
)::TracerAuxiliaries
    return TracerAuxiliaries(
        [zeros(0, 0, 0) for field in fieldnames(TracerAuxiliaries)]...,
    )
end

function TracerAuxiliaries(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    tracer_setup::Val{:TracerOn},
    variables::Variables,
)::TracerAuxiliaries
    (; nxx, nyy, nzz, i0, i1, j0, j1) = domain
    (; x, y, zc) = grid
    (; lref) = constants
    (; rhobar) = atmosphere
    (; lref) = constants
    (; background_tracer) = namelists.tracer

    backgroundtracer = zeros(nxx, nyy, nzz)
    @ivy for k in 1:nzz, j in j0:j1, i in i0:i1
        backgroundtracer[i, j, k] =
            background_tracer(x[i] * lref, y[j] * lref, zc[i, j, k] * lref)
    end
    set_zonal_boundaries_of_field!(backgroundtracer, namelists, domain)
    set_meridional_boundaries_of_field!(backgroundtracer, namelists, domain)

    backgroundtracer .*= rhobar

    return TracerAuxiliaries(backgroundtracer)
end
