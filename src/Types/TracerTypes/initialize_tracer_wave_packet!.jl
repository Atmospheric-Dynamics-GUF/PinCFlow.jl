"""
```julia
initialize_tracer_wave_packet!(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
    alphatracer::AbstractFloat,
    chi::AbstractArray{<:AbstractFloat, 3},
    testcase::AbstractTestCase,
)
```

Return for non-wave-packet test cases.

```julia
initialize_tracer_wave_packet!(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
    alphatracer::AbstractFloat,
    chi::AbstractArray{<:AbstractFloat, 3},
    testcase::WavePacket,
)
```

Initialize the fluctuations of a non-dimensional tracer with a wave-packet inferred from the density fluctuations.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `atmosphere`: Atmospheric background fields.

  - `grid`: Collection of parameters and fields that describe the grid.

  - `variables`: Container for arrays needed for the prediction of the prognostic variables.

  - `alphatracer`: Ratio between the tracer fluctuations and the buoyancy fluctuations divided by the squared buoyancy frequency.

  - `chi`: Non-dimensional tracer.

  - `testcase`: Test case on which the current simulation is based.
"""
function initialize_tracer_wave_packet! end

function initialize_tracer_wave_packet!(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
    alphatracer::AbstractFloat,
    chi::AbstractArray{<:AbstractFloat, 3},
    testcase::AbstractTestCase,
)
    return
end

function initialize_tracer_wave_packet!(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
    alphatracer::AbstractFloat,
    chi::AbstractArray{<:AbstractFloat, 3},
    testcase::WavePacket,
)
    (; rhop) = variables.predictands
    (; bvsstrattfc, rhostrattfc) = atmosphere
    (; fr2) = constants
    (; nxx, nyy, nzz) = domain

    bprime = (rhostrattfc ./ (rhop .+ rhostrattfc) .- 1.0) ./ fr2

    chiprime = zeros(nxx, nyy, nzz)

    chiprime = alphatracer ./ bvsstrattfc .* bprime

    chi .= chi .+ chiprime

    return
end
