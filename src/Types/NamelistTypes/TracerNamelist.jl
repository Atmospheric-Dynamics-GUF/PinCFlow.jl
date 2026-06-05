"""
```julia
TracerNamelist
```

Namelist for the inclusion of a tracer and the calculation of the leading-order gravity-wave impact.

```julia
TracerNamelist(;
    tracer_setup::Symbol = :NoTracer,
    leading_order_impact::Bool = false,
    initial_chi::Function = (x, y, z) -> 0.0,
    relaxed_chi::Function = (x, y, z, t, dt) -> 0.0,
    apply_lhs_sponge_to_tracer::Bool = true,
)::TracerNamelist
```

Construct a `TracerNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `tracer_setup::Symbol`: General tracer configuration.

  - `leading_order_impact::Bool`: Flag to include the leading-order impact of gravity waves when parameterizing waves with the WKB model.

  - `initial_chi::FunctionWrapper{Float64, NTuple{3, Float64}}`: Function used to initialize the tracer.

  - `relaxed_chi::FunctionWrapper{Float64, NTuple{5, Float64}}`: Function used to compute the background tracer.

  - `apply_lhs_sponge_to_tracer::Bool`: Flag to relax the tracer fields to `relaxed_chi`.
"""
struct TracerNamelist
    tracer_setup::Symbol
    leading_order_impact::Bool
    initial_chi::FunctionWrapper{Float64, NTuple{3, Float64}}
    relaxed_chi::FunctionWrapper{Float64, NTuple{5, Float64}}
    apply_lhs_sponge_to_tracer::Bool
end

function TracerNamelist(;
    tracer_setup::Symbol = :NoTracer,
    leading_order_impact::Bool = false,
    initial_chi::Function = (x, y, z) -> 0.0,
    relaxed_chi::Function = (x, y, z, t, dt) -> 0.0,
    apply_lhs_sponge_to_tracer::Bool = true,
)::TracerNamelist
    return TracerNamelist(
        tracer_setup,
        leading_order_impact,
        initial_chi,
        relaxed_chi,
        apply_lhs_sponge_to_tracer,
    )
end
