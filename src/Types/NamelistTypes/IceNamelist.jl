"""
```julia
IceNamelist{A <: AbstractIce}
```

Namelist for the inclusion of ice physics.

```julia
IceNamelist(; icesetup::AbstractIce = NoIce())
```

Construct an `IceNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `icesetup::A`: General ice-physics configuration.
"""
struct IceNamelist{A <: AbstractIce}
    icesetup::A
end

function IceNamelist(; icesetup::AbstractIce = NoIce())
    return IceNamelist(icesetup)
end
