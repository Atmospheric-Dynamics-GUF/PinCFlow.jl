"""
```julia
TracerNamelist{A <: AbstractTracer, B <: Bool}
```

Namelist for the inclusion of a tracer and the calculation of the leading-order gravity-wave impact.

```julia
TracerNamelist(;
    tracersetup::AbstractTracer = NoTracer(),
    leading_order_impact = true,
)::TracerNamelist
```

Construct a `TracerNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `tracersetup::A`: General tracer configuration.

  - `leading_order_impact::B`: Flag to include the leading-order impact of gravity waves when parameterizing waves with the WKB model.
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
