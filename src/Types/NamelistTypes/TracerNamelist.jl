"""
```julia
TracerNamelist
```

Namelist for the inclusion of a tracer and the calculation of the leading-order gravity-wave impact.

```julia
TracerNamelist(;
    tracer_setup::Symbol = :NoTracer,
    leading_order_impact::Bool = false,
    initial_tracer::Function = (x, y, z) -> 0.0,
    background_tracer::Function = (x, y, z) -> 0.0,
    apply_sponge_to_tracer::Bool = true,
)::TracerNamelist
```

Construct a `TracerNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `tracer_setup::Symbol`: General tracer configuration.

  - `leading_order_impact::Bool`: Flag to include the leading-order impact of gravity waves when parameterizing waves with the WKB model.

  - `initial_tracer::FunctionWrapper{Float64, NTuple{3, Float64}}`: Function used to initialize the tracer.

  - `background_tracer::FunctionWrapper{Float64, NTuple{3, Float64}}`: Function used to compute the background tracer.

  - `apply_sponge_to_tracer::Bool`: Flag to relax the tracer fields to `backgroundtracer`.
"""
struct TracerNamelist
    tracer_setup::Symbol
    leading_order_impact::Bool
    initial_tracer::FunctionWrapper{Float64, NTuple{3, Float64}}
    background_tracer::FunctionWrapper{Float64, NTuple{3, Float64}}
    apply_sponge_to_tracer::Bool
end

function TracerNamelist(;
    tracer_setup::Symbol = :NoTracer,
    leading_order_impact::Bool = false,
    initial_tracer::Function = (x, y, z) -> 0.0,
    background_tracer::Function = (x, y, z) -> 0.0,
    apply_sponge_to_tracer::Bool = true,
)::TracerNamelist
    return TracerNamelist(tracer_setup, leading_order_impact, initial_tracer, background_tracer, apply_sponge_to_tracer)
end
