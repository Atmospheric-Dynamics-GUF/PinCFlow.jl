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
struct IceNamelist{A <: AbstractIce, B <: AbstractFloat, C <: Integer, D <: AbstractCloudCover}
     icesetup::A
     dt_ice::B
     nscx :: C
     nscy :: C
     nscz :: C
     cloudcover :: D
     parameterized_nucleation :: Bool
     parameterized_sgs_q :: Bool
end

function IceNamelist(; icesetup::AbstractIce, dt_ice = 1.0, nscx = 1, nscy = 1, nscz = 1, cloudcover= CloudCoverOff(), parameterized_nucleation = false, parameterized_sgs_q = false)
    return IceNamelist(icesetup, dt_ice, nscx, nscy, nscz, cloudcover, parameterized_nucleation, parameterized_sgs_q)
end

