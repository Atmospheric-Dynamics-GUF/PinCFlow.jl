"""
```julia
TracerNamelist
```

Namelist for the inclusion of a tracer and the calculation of the leading-order gravity-wave impact.

```julia
TracerNamelist(;
    tracer_setup::Symbol = :no_tracer,
    leading_order_impact::Bool = false,
    initial_tracer::Function = (x, y, z) -> 0.0,
)::TracerNamelist
```

Construct a `TracerNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `tracer_setup::Symbol`: General tracer configuration.

  - `leading_order_impact::Bool`: Flag to include the leading-order impact of gravity waves when parameterizing waves with the WKB model.

  - `initial_tracer::FunctionWrapper{Float64, NTuple{3, Float64}}`: Function used to initialize the tracer.
"""
struct TracerNamelist
    tracer_setup::Symbol
    leading_order_impact::Bool
    initial_tracer::FunctionWrapper{Float64, NTuple{3, Float64}}
end

function TracerNamelist(;
    tracer_setup::Symbol = :no_tracer,
    leading_order_impact::Bool = false,
    initial_tracer::Function = (x, y, z) -> 0.0,
)::TracerNamelist
    return TracerNamelist(tracer_setup, leading_order_impact, initial_tracer)
end
