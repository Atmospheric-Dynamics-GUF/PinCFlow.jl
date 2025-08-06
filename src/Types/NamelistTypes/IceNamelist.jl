"""
```julia
IceNamelist{A <: AbstractIce}
```

```julia
IceNamelist(; icesetup::AbstractIce = NoIce())
```
"""
struct IceNamelist{A <: AbstractIce}
    icesetup::A
end

function IceNamelist(; icesetup::AbstractIce = NoIce())
    return IceNamelist(icesetup)
end
