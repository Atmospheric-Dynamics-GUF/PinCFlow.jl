"""
```julia
TracerForcings{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
}
```

```julia
TracerForcings(namelists::Namelists, domain::Domain)::TracerForcings
```
"""
struct TracerForcings{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
}
    chiu0::A
    chiv0::A
    chiw0::A
    chiu1::B
    chiv1::B
    chiw1::B
    chiturbu::C
    chiturbv::C
    chiturbw::C
end

function TracerForcings(namelists::Namelists, domain::Domain)::TracerForcings
    (; tracersetup) = namelists.tracer

    return
end
