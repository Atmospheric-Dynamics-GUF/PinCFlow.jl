"""
```julia
IceNamelist{A <: AbstractIce}
```
"""
struct IceNamelist{A <: AbstractIce}
    icesetup::A
end

"""
```julia
IceNamelist(; icesetup = NoIce())
```
"""
function IceNamelist(; icesetup = NoIce())
    return IceNamelist(icesetup)
end
