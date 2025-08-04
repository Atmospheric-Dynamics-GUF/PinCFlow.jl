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
TracerNamelist(; tracersetup::AbstractTracer = NoTracer())
```
"""
function TracerNamelist(; tracersetup::AbstractTracer = NoTracer())
    return TracerNamelist(tracersetup)
end
