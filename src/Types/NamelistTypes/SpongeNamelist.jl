"""
```julia
SpongeNamelist
```

Namelist for parameters describing the sponge.

```julia
SpongeNamelist(;
    damp_horizontal_wind_on_rhs::Bool = false,
    relax_to_mean::Bool = false,
    lhs_sponge::Function = (x, y, z, t, dt) -> 0.0,
    rhs_sponge::Function = (x, y, z, t, dt) -> 0.0,
    relaxed_u::Function = (x, y, z, t, dt) -> 0.0,
    relaxed_v::Function = (x, y, z, t, dt) -> 0.0,
    relaxed_w::Function = (x, y, z, t, dt) -> 0.0,
)::SpongeNamelist
```

Construct a `SpongeNamelist` instance with the given keyword arguments as properties, converting them to meet the type constraints.

# Fields/Keywords

  - `damp_horizontal_wind_on_rhs::Bool`: Switch for applying the RHS sponge to the horizontal wind.

  - `relax_to_mean::Bool`: Switch for relaxing the wind towards its averages on the terrain-following surfaces. If set to `false`, the relaxation wind is computed with `relaxed_u`, `relaxed_v` and `relaxed_w`.

  - `lhs_sponge::FunctionWrapper{Float64, NTuple{5, Float64}}`: Function used to compute the Rayleigh-damping coefficient of the LHS sponge.

  - `rhs_sponge::FunctionWrapper{Float64, NTuple{5, Float64}}`: Function used to compute the Rayleigh-damping coefficient of the RHS sponge.

  - `relaxed_u::FunctionWrapper{Float64, NTuple{5, Float64}}`: Function used to compute the zonal relaxation wind if `relax_to_mean` is set to `false`.

  - `relaxed_v::FunctionWrapper{Float64, NTuple{5, Float64}}`: Function used to compute the meridional relaxation wind if `relax_to_mean` is set to `false`.

  - `relaxed_w::FunctionWrapper{Float64, NTuple{5, Float64}}`: Function used to compute the vertical relaxation wind if `relax_to_mean` is set to `false`.
"""
struct SpongeNamelist
    damp_horizontal_wind_on_rhs::Bool
    relax_to_mean::Bool
    lhs_sponge::FunctionWrapper{Float64, NTuple{5, Float64}}
    rhs_sponge::FunctionWrapper{Float64, NTuple{5, Float64}}
    relaxed_u::FunctionWrapper{Float64, NTuple{5, Float64}}
    relaxed_v::FunctionWrapper{Float64, NTuple{5, Float64}}
    relaxed_w::FunctionWrapper{Float64, NTuple{5, Float64}}
end

function SpongeNamelist(;
    damp_horizontal_wind_on_rhs::Bool = false,
    relax_to_mean::Bool = false,
    lhs_sponge::Function = (x, y, z, t, dt) -> 0.0,
    rhs_sponge::Function = (x, y, z, t, dt) -> 0.0,
    relaxed_u::Function = (x, y, z, t, dt) -> 0.0,
    relaxed_v::Function = (x, y, z, t, dt) -> 0.0,
    relaxed_w::Function = (x, y, z, t, dt) -> 0.0,
)::SpongeNamelist
    return SpongeNamelist(
        damp_horizontal_wind_on_rhs,
        relax_to_mean,
        lhs_sponge,
        rhs_sponge,
        relaxed_u,
        relaxed_v,
        relaxed_w,
    )
end
