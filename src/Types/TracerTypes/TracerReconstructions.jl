"""
```julia
TracerReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
```
"""
struct TracerReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
    chitilde::A
end

"""
```julia
TracerReconstructions(namelists::Namelists, domain::Domain)
```
"""
function TracerReconstructions(namelists::Namelists, domain::Domain)
    (; tracersetup) = namelists.tracer

    return TracerReconstructions(domain, tracersetup)
end

"""
```julia
TracerReconstructions(domain::Domain, tracersetup::NoTracer)
```
"""
function TracerReconstructions(domain::Domain, tracersetup::NoTracer)
    chitilde = zeros(0, 0, 0, 0, 0)

    return TracerReconstructions(chitilde)
end

"""
```julia
TracerReconstructions(domain::Domain, tracersetup::AbstractTracer)
```
"""
function TracerReconstructions(domain::Domain, tracersetup::AbstractTracer)
    (; nxx, nyy, nzz) = domain

    chitilde = zeros(nxx, nyy, nzz, 3, 2)

    return TracerReconstructions(chitilde)
end
