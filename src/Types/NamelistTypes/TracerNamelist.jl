"""
```julia
TracerNamelist{A <: AbstractTracer, B <: Bool, C <: Function}
```

Namelist for the inclusion of a tracer and the calculation of the leading-order gravity-wave impact.

```julia
TracerNamelist(;
    tracer_setup::AbstractTracer = NoTracer(),
    leading_order_impact::Bool = false,
    initial_tracer::Function = (x, y, z) -> 0.0,
)::TracerNamelist
```

Construct a `TracerNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `tracer_setup::A`: General tracer configuration.

  - `leading_order_impact::B`: Flag to include the leading-order impact of gravity waves when parameterizing waves with the WKB model.

  - `initial_tracer::C`: Function used to initialize the tracer.
"""
struct TracerNamelist{A <: AbstractTracer, B <: Bool, C <: Function}
    tracer_setup::A
    leading_order_impact::B
    initial_tracer::C
end

function TracerNamelist(;
    tracer_setup::AbstractTracer = NoTracer(),
    leading_order_impact::Bool = false,
    initial_tracer::Function = (x, y, z) -> 0.0,
)::TracerNamelist
    return TracerNamelist(tracer_setup, leading_order_impact, initial_tracer)
end
