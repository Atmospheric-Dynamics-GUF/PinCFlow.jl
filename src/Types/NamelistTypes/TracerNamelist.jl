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
struct TracerNamelist{A <: AbstractTracer, B <: Bool}
    tracersetup::A
    leading_order_impact::B
end

function TracerNamelist(;

    tracersetup::AbstractTracer = NoTracer(),
    leading_order_impact = true,
)::TracerNamelist
    return TracerNamelist(tracersetup, leading_order_impact)
end
