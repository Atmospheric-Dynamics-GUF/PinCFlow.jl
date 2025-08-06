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
)
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
)
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
)
```

Construct a `TracerPredictands` instance with an initialized non-dimensional tracer.

The initialization consists of three steps. First, the tracer is set to ``\\chi = z``, then fluctuations in the shape of a wave packet are added if the test case calls for it (see `initialize_tracer_wave_packet!`). Finally, the density is absorbed into the tracer, i.e. the substitution ``\\chi \\rightarrow \\rho \\chi`` is performed.

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

# See also

  - [`PinCFlow.Types.TracerTypes.initialize_tracer_wave_packet!`](@ref)
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
)
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
)
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
)
    (; nxx, nyy, nzz) = domain
    (; ztfc) = grid
    (; rhostrattfc) = atmosphere
    (; rho) = variables.predictands
    (; testcase) = namelists.setting

    alphatracer = 1.0
    chi = zeros(nxx, nyy, nzz)
    chi .= alphatracer .* ztfc

    initialize_tracer_wave_packet!(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        variables,
        alphatracer,
        chi,
        testcase,
    )

    chi .= chi .* (rho .+ rhostrattfc)

    return TracerPredictands(chi)
end
