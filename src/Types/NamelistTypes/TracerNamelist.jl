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
    background_tracer::Function = (x, y, z) -> 0.0,
)::TracerNamelist
```

Construct a `TracerNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `tracer_setup::A`: General tracer configuration.

  - `leading_order_impact::B`: Flag to include the leading-order impact of gravity waves when parameterizing waves with the WKB model.

  - `initial_tracer::C`: Function used to initialize the tracer.

  - `background_tracer::C`: Function used to set the background tracer distribution.
"""
struct TracerNamelist{A <: AbstractTracer, B <: Bool, C <: Function, D <: Function}
    tracer_setup::A
    leading_order_impact::B
    next_order_impact::B
    turbulence_impact::B
    initial_tracer::C
    background_tracer::D
end

function TracerNamelist(;
    tracer_setup::AbstractTracer = NoTracer(),
    leading_order_impact::Bool = true,
    next_order_impact::Bool = true,
    turbulence_impact::Bool = true,
    initial_tracer::Function = (x, y, z) -> 0.0,
    background_tracer::Function = (x, y, z) -> 0.0,
)::TracerNamelist
    return TracerNamelist(
        tracer_setup,
        leading_order_impact,
        next_order_impact,
        turbulence_impact,
        initial_tracer,
        background_tracer,
    )
end
