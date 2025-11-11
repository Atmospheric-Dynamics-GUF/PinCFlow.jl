"""
```julia
TracerNamelist{A <: AbstractTracer, B <: Bool, C <: AbstractFloat}
```

Namelist for the inclusion of a tracer and the calculation of the leading-order gravity-wave impact.

```julia
TracerNamelist(;
    tracersetup::AbstractTracer = NoTracer(),
    leading_order_impact::Bool = false,
    alphatracer::AbstractFloat = 1.0E+0,
)::TracerNamelist
```

Construct a `TracerNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `tracersetup::A`: General tracer configuration.

  - `leading_order_impact::B`: Flag to include the leading-order impact of gravity waves when parameterizing waves with the WKB model.

  - `alphatracer::C`: Scaling parameter for the initial tracer distribution.
"""
struct TracerNamelist{A <: AbstractTracer, B <: Bool, C <: AbstractFloat}
    tracersetup::A
    leading_order_impact::B
    alphatracer::C
end

function TracerNamelist(;
    tracersetup::AbstractTracer = NoTracer(),
    leading_order_impact::Bool = false,
    alphatracer::AbstractFloat = 1.0E+0,
)::TracerNamelist
    return TracerNamelist(tracersetup, leading_order_impact, alphatracer)
end
