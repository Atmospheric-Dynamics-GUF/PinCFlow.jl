"""
```julia
DiscretizationNamelist
```

Namelist for parameters describing the discretization.

```julia
DiscretizationNamelist(;
    cfl_number::Real = 5.0E-1,
    wkb_cfl_number::Real = 5.0E-1,
    dtmin::Real = 1.0E-6,
    dtmax::Real = 1.0E+3,
    adaptive_time_step::Bool = true,
    limiter_type::Symbol = :monotone_centered_variant,
)::DiscretizationNamelist
```

Construct a `DiscretizationNamelist` instance with the given keyword arguments as properties, converting them to meet the type constraints.

# Fields/Keywords

  - `cfl_number::Float64`: Number used for the CFL condition in the time step computation.

  - `wkb_cfl_number::Float64`: Number used for the WKB-CFL condition in the time step computation.

  - `dtmin::Float64`: Minimum time step allowed for the integration.

  - `dtmax::Float64`: Maximum time step allowed for the integration.

  - `adaptive_time_step::Bool`: Switch for using stability criteria to determine the time step. If set to `false`, `dtmax` is used as a fixed time step.

  - `limiter_type::Symbol`: Flux limiter used by the MUSCL scheme.
"""
struct DiscretizationNamelist
    cfl_number::Float64
    wkb_cfl_number::Float64
    dtmin::Float64
    dtmax::Float64
    adaptive_time_step::Bool
    limiter_type::Symbol
end

function DiscretizationNamelist(;
    cfl_number::Real = 5.0E-1,
    wkb_cfl_number::Real = 5.0E-1,
    dtmin::Real = 1.0E-6,
    dtmax::Real = 1.0E+3,
    adaptive_time_step::Bool = true,
    limiter_type::Symbol = :monotone_centered_variant,
)::DiscretizationNamelist
    return DiscretizationNamelist(
        Float64(cfl_number),
        Float64(wkb_cfl_number),
        Float64(dtmin),
        Float64(dtmax),
        adaptive_time_step,
        limiter_type,
    )
end
