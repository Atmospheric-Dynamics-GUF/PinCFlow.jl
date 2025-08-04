"""
```julia
TracerFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
```
"""
struct TracerFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
    phichi::A
end

"""
```julia
TracerFluxes(namelists::Namelists, domain::Domain)
```
"""
function TracerFluxes(namelists::Namelists, domain::Domain)
    (; tracersetup) = namelists.tracer

    return TracerFluxes(domain, tracersetup)
end

"""
```julia
TracerFluxes(domain::Domain, tracersetup::NoTracer)
```
"""
function TracerFluxes(domain::Domain, tracersetup::NoTracer)
    phichi = zeros(0, 0, 0, 0)

    return TracerFluxes(phichi)
end

"""
```julia
TracerFluxes(domain::Domain, tracersetup::AbstractTracer)
```
"""
function TracerFluxes(domain::Domain, tracersetup::AbstractTracer)
    (; nxx, nyy, nzz) = domain

    phichi = zeros(nxx, nyy, nzz, 3)

    return TracerFluxes(phichi)
end
