"""
```julia
IceNamelist{A <: AbstractIce}
```

```julia
IceNamelist(; icesetup::AbstractIce = NoIce())
```
"""
struct IceNamelist{A <: AbstractIce, B <: AbstractFloat}
     icesetup::A
     dt_ice::B
end

function IceNamelist(; icesetup::AbstractIce = NoIce(), dt_ice = 1.0)
    return IceNamelist(icesetup, dt_ice)
end

