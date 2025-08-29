"""
```julia
IceNamelist{A <: AbstractIce}
```

```julia
IceNamelist(; icesetup::AbstractIce = NoIce())
```
"""
struct IceNamelist{A <: AbstractIce, B <: AbstractFloat, C <: Integer}
     icesetup::A
     dt_ice::B
     compute_cloudcover :: C
     nscx :: C
     nscy :: C
     nscz :: C
end

function IceNamelist(; icesetup::AbstractIce = NoIce(), dt_ice = 1.0, compute_cloudcover = 0, nscx = 1, nscy = 1, nscz = 1)
    return IceNamelist(icesetup, dt_ice, compute_cloudcover, nscx, nscy, nscz)
end

