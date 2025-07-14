"""
```julia
DiscretizationNamelist{A <: AbstractFloat, B <: Bool, C <: AbstractLimiter}
```

Namelist for the discretization (see constructor for parameter descriptions).
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

"""
```julia
DiscretizationNamelist(;
    cfl = 5.0E-1,
    cfl_wave = 5.0E-1,
    dtmin_dim = 1.0E-6,
    dtmax_dim = 1.0E+3,
    adaptive_time_step = true,
    limitertype = MCVariant(),
)
```

Construct a DiscretizationNamelist object, which holds parameters for the time discretization.

# Arguments

  - `cfl::Float = 0.5`: Courant-Friedrichs-Lewy number used in `compute_time_step`.
  - `cfl_wave::Float = 0.5`: Courant-Friedrichs-Lewy number for WKB mode.
  - `dtmin_dim::Float = 1.0E-6`: minimum time step the integration routine is allowed to take.
  - `dtmax_dim::Float = 1.0E+3`: maximum time step the integration routine is allowed to take.
  - `adaptive_time_step::Bool = true`: whether to use adaptive time step. If false, a constant time step based on `dtmax_dim` is used.
  - `limitertype::AbstractLimiter = MCVariant()`: type of limiter used in MUSCL scheme. Currently the only option is `MCVariant()`.
"""
function DiscretizationNamelist(;
    cfl = 5.0E-1,
    cfl_wave = 5.0E-1,
    dtmin_dim = 1.0E-6,
    dtmax_dim = 1.0E+3,
    adaptive_time_step = true,
    limitertype = MCVariant(),
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
