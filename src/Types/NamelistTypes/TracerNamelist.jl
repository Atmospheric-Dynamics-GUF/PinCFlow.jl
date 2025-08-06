"""
```julia
TracerNamelist{A <: AbstractTracer}
```

```julia
TracerNamelist(; tracersetup::AbstractTracer = NoTracer())
```
"""
struct TracerNamelist{A <: AbstractTracer}
    tracersetup::A
end

function TracerNamelist(; tracersetup::AbstractTracer = NoTracer())
    return TracerNamelist(tracersetup)
end
