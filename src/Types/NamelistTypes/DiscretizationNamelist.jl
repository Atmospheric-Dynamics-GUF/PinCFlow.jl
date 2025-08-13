"""
```julia
DiscretizationNamelist{A <: AbstractFloat, B <: Bool, C <: AbstractLimiter}
```

Namelist for parameters describing the discretization.

```julia
DiscretizationNamelist(;
    cfl::AbstractFloat = 5.0E-1,
    cfl_wave::AbstractFloat = 5.0E-1,
    dtmin_dim::AbstractFloat = 1.0E-6,
    dtmax_dim::AbstractFloat = 1.0E+3,
    adaptive_time_step::Bool = true,
    limitertype::AbstractLimiter = MCVariant(),
)
```

Construct a `DiscretizationNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `cfl::A`: Number used for the CFL condition in the time step computation.

  - `cfl_wave::A`: Number used for the WKB-CFL condition in the time step computation.

  - `dtmin_dim::A`: Minimum time step allowed for the integration.

  - `dtmax_dim::A`: Maximum time step allowed for the integration.

  - `adaptive_time_step::B`: Switch for using stability criteria to determine the time step. If set to `false`, `dtmax_dim` is used as a fixed time step.

  - `limitertype::C`: Flux limiter used by the MUSCL scheme.
"""
struct DiscretizationNamelist{
    A <: AbstractFloat,
    B <: Bool,
    C <: AbstractLimiter,
}
    cfl::A
    cfl_wave::A
    dtmin_dim::A
    dtmax_dim::A
    adaptive_time_step::B
    limitertype::C
end

function DiscretizationNamelist(;
    cfl::AbstractFloat = 5.0E-1,
    cfl_wave::AbstractFloat = 5.0E-1,
    dtmin_dim::AbstractFloat = 1.0E-6,
    dtmax_dim::AbstractFloat = 1.0E+3,
    adaptive_time_step::Bool = true,
    limitertype::AbstractLimiter = MCVariant(),
)
    return DiscretizationNamelist(
        cfl,
        cfl_wave,
        dtmin_dim,
        dtmax_dim,
        adaptive_time_step,
        limitertype,
    )
end
