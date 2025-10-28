"""
```julia
DiscretizationNamelist{
    A <: AbstractFloat,
    B <: AbstractFloat,
    C <: AbstractFloat,
    D <: AbstractFloat,
    E <: Bool,
    F <: AbstractLimiter,
}
```

Namelist for parameters describing the discretization.

```julia
DiscretizationNamelist(;
    cfl_number::AbstractFloat = 5.0E-1,
    wkb_cfl_number::AbstractFloat = 5.0E-1,
    dtmin::AbstractFloat = 1.0E-6,
    dtmax::AbstractFloat = 1.0E+3,
    adaptive_time_step::Bool = true,
    limiter_type::AbstractLimiter = MCVariant(),
)::DiscretizationNamelist
```

Construct a `DiscretizationNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `cfl_number::A`: Number used for the CFL condition in the time step computation.

  - `wkb_cfl_number::B`: Number used for the WKB-CFL condition in the time step computation.

  - `dtmin::C`: Minimum time step allowed for the integration.

  - `dtmax::D`: Maximum time step allowed for the integration.

  - `adaptive_time_step::E`: Switch for using stability criteria to determine the time step. If set to `false`, `dtmax` is used as a fixed time step.

  - `limiter_type::F`: Flux limiter used by the MUSCL scheme.
"""
struct DiscretizationNamelist{
    A <: AbstractFloat,
    B <: AbstractFloat,
    C <: AbstractFloat,
    D <: AbstractFloat,
    E <: Bool,
    F <: AbstractLimiter,
}
    cfl_number::A
    wkb_cfl_number::B
    dtmin::C
    dtmax::D
    adaptive_time_step::E
    limiter_type::F
end

function DiscretizationNamelist(;
    cfl_number::AbstractFloat = 5.0E-1,
    wkb_cfl_number::AbstractFloat = 5.0E-1,
    dtmin::AbstractFloat = 1.0E-6,
    dtmax::AbstractFloat = 1.0E+3,
    adaptive_time_step::Bool = true,
    limiter_type::AbstractLimiter = MCVariant(),
)::DiscretizationNamelist
    return DiscretizationNamelist(
        cfl_number,
        wkb_cfl_number,
        dtmin,
        dtmax,
        adaptive_time_step,
        limiter_type,
    )
end
