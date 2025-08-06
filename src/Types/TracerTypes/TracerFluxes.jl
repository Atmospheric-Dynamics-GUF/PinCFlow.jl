"""
```julia
TracerFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
```

```julia
TracerFluxes(namelists::Namelists, domain::Domain)
```

```julia
TracerFluxes(domain::Domain, tracersetup::NoTracer)
```

```julia
TracerFluxes(domain::Domain, tracersetup::AbstractTracer)
```
"""
struct TracerFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
    phichi::A
end

function TracerFluxes(namelists::Namelists, domain::Domain)
    (; tracersetup) = namelists.tracer

    return TracerFluxes(domain, tracersetup)
end

function TracerFluxes(domain::Domain, tracersetup::NoTracer)
    phichi = zeros(0, 0, 0, 0)

    return TracerFluxes(phichi)
end

function TracerFluxes(domain::Domain, tracersetup::AbstractTracer)
    (; nxx, nyy, nzz) = domain

    phichi = zeros(nxx, nyy, nzz, 3)

    return TracerFluxes(phichi)
end
