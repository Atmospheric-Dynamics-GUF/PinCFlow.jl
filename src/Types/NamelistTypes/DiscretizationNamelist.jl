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
    cfl::AbstractFloat = 5.0E-1,
    cfl_wave::AbstractFloat = 5.0E-1,
    dtmin_dim::AbstractFloat = 1.0E-6,
    dtmax_dim::AbstractFloat = 1.0E+3,
    adaptive_time_step::Bool = true,
    limitertype::AbstractLimiter = MCVariant(),
)
```

Construct a DiscretizationNamelist object, which holds parameters for the time discretization.

# Arguments

  - `cfl`: Courant-Friedrichs-Lewy number used in `compute_time_step`.
  - `cfl_wave`: Courant-Friedrichs-Lewy number for WKB mode.
  - `dtmin_dim`: minimum time step the integration routine is allowed to take.
  - `dtmax_dim`: maximum time step the integration routine is allowed to take.
  - `adaptive_time_step`: whether to use adaptive time step. If false, a constant time step based on `dtmax_dim` is used.
  - `limitertype`: type of limiter used in MUSCL scheme. Currently the only option is `MCVariant()`.

# Returns

  - `::DiscretizationNamelist`: `DiscretizationNamelist` instance.
"""
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
