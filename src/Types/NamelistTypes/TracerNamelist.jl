"""
```julia
TracerNamelist{A <: AbstractTracer, B <: Bool, C <: AbstractFloat}
```

Namelist for the inclusion of a tracer and the calculation of the leading-order gravity-wave impact.

```julia
TracerNamelist(;
    tracer_setup::AbstractTracer = NoTracer(),
    leading_order_impact::Bool = false,
    alphatracer::AbstractFloat = 1.0E+0,
)::TracerNamelist
```

Construct a `TracerNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `tracer_setup::A`: General tracer configuration.

  - `leading_order_impact::B`: Flag to include the leading-order impact of gravity waves when parameterizing waves with the WKB model.

  - `alphatracer::C`: Scaling parameter for the initial tracer distribution.
"""
struct TracerNamelist{A <: AbstractTracer, B <: Bool, C <: AbstractFloat}
    tracer_setup::A
    leading_order_impact::B
    alphatracer::C
end

function TracerNamelist(;
    tracer_setup::AbstractTracer = NoTracer(),
    leading_order_impact::Bool = false,
    alphatracer::AbstractFloat = 1.0E+0,
)::TracerNamelist
    return TracerNamelist(tracer_setup, leading_order_impact, alphatracer)
end
