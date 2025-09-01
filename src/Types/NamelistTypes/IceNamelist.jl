"""
```julia
IceNamelist{A <: AbstractIce}
```

Namelist for the inclusion of ice physics.

```julia
IceNamelist(; icesetup::AbstractIce = NoIce())::IceNamelist
```

Construct an `IceNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `icesetup::A`: General ice-physics configuration.
"""
struct IceNamelist{A <: AbstractIce, B <: AbstractFloat, C <: Integer}
     icesetup::A
     dt_ice::B
     compute_cloudcover :: C
     nscx :: C
     nscy :: C
     nscz :: C
end

function IceNamelist(; icesetup::AbstractIce, dt_ice = 1.0, compute_cloudcover = 0, nscx = 1, nscy = 1, nscz = 1)
    return IceNamelist(icesetup, dt_ice, compute_cloudcover, nscx, nscy, nscz)
end

