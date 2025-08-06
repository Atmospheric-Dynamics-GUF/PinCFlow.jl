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
