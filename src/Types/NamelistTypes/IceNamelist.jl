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
IceNamelist(; icesetup::AbstractIce = NoIce())
```
"""
function IceNamelist(; icesetup::AbstractIce = NoIce())
    return IceNamelist(icesetup)
end
