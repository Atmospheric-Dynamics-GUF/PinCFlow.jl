"""
```julia
TracerAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```

Initial states of the tracers.

```julia
TracerAuxiliaries(tracerpredictands::TracerPredictands)::TracerAuxiliaries
```

Construct a `TracerAuxiliaries` instance by copying the arrays in `tracerpredictands`.

# Fields

  - `backgroundtracer::A`: Initial state of a non-dimensional tracer.

# Arguments

  - `tracerpredictands`: Tracers.
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
    tracer_setup::NoTracer,
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
    tracer_setup::TracerOn,
    variables::Variables,
)::TracerAuxiliaries
    (; nxx, nyy, nzz, i0, i1, j0, j1) = domain
    (; x, y, zc) = grid
    (; lref) = constants
    (; rhobar) = atmosphere
    (; rho) = variables.predictands
    (; lref) = constants
    (; background_tracer) = namelists.tracer

    backgroundtracer = zeros(nxx, nyy, nzz)
    @ivy for k in 1:nzz, j in j0:j1, i in i0:i1
        backgroundtracer[i, j, k] =
            background_tracer(x[i] * lref, y[j] * lref, zc[i, j, k] * lref)
    end
    set_zonal_boundaries_of_field!(backgroundtracer, namelists, domain)
    set_meridional_boundaries_of_field!(backgroundtracer, namelists, domain)

    backgroundtracer .*= rho .+ rhobar

    return TracerAuxiliaries(backgroundtracer)
end
