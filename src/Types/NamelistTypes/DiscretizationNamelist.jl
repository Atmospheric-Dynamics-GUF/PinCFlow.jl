"""
```julia
DiscretizationNamelist{A <: AbstractFloat, B <: Bool, C <: AbstractLimiter}
```

Namelist for parameters describing the discretization.

```julia
DiscretizationNamelist(;
    cfl_number::AbstractFloat = 5.0E-1,
    wkb_cfl_number::AbstractFloat = 5.0E-1,
    dtmin::AbstractFloat = 1.0E-6,
    dtmax::AbstractFloat = 1.0E+3,
    adaptive_time_step::Bool = true,
    limitertype::AbstractLimiter = MCVariant(),
)::DiscretizationNamelist
```

Construct a `DiscretizationNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `cfl_number::A`: Number used for the CFL condition in the time step computation.

  - `wkb_cfl_number::A`: Number used for the WKB-CFL condition in the time step computation.

  - `dtmin::A`: Minimum time step allowed for the integration.

  - `dtmax::A`: Maximum time step allowed for the integration.

  - `adaptive_time_step::B`: Switch for using stability criteria to determine the time step. If set to `false`, `dtmax` is used as a fixed time step.

  - `limitertype::C`: Flux limiter used by the MUSCL scheme.
"""
struct DiscretizationNamelist{
    A <: AbstractFloat,
    B <: Bool,
    C <: AbstractLimiter,
}
    cfl_number::A
    wkb_cfl_number::A
    dtmin::A
    dtmax::A
    adaptive_time_step::B
    limitertype::C
end

function DiscretizationNamelist(;
    cfl_number::AbstractFloat = 5.0E-1,
    wkb_cfl_number::AbstractFloat = 5.0E-1,
    dtmin::AbstractFloat = 1.0E-6,
    dtmax::AbstractFloat = 1.0E+3,
    adaptive_time_step::Bool = true,
    limitertype::AbstractLimiter = MCVariant(),
)::DiscretizationNamelist
    return DiscretizationNamelist(
        cfl_number,
        wkb_cfl_number,
        dtmin,
        dtmax,
        adaptive_time_step,
        limitertype,
    )
end
