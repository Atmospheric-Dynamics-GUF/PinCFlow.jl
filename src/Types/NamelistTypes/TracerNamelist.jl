"""
```julia
TracerNamelist{A <: AbstractTracer}
```
"""
struct TracerNamelist{A <: AbstractTracer}
    tracersetup::A
end

"""
```julia
TracerNamelist(; tracersetup = NoTracer())
```
"""
function TracerNamelist(; tracersetup = NoTracer())
    return TracerNamelist(tracersetup)
end
