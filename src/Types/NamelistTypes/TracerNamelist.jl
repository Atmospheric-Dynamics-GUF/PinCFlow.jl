"""
```julia
TracerNamelist{A <: AbstractTracer}
```

Namelist for the inclusion of a tracer.

```julia
TracerNamelist(;
    tracersetup::AbstractTracer = NoTracer(),
)::TracerNamelist
```

Construct a `TracerNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `tracersetup::A`: General tracer configuration.
"""
struct TracerNamelist{A <: AbstractTracer}
    tracersetup::A
end

function TracerNamelist(;
    tracersetup::AbstractTracer = NoTracer(),
)::TracerNamelist
    return TracerNamelist(tracersetup)
end
