"""
```julia
SpongeNamelist{
    A <: Bool,
    B <: Function,
    C <: Function,
    D <: Function,
    E <: Function,
    F <: Function,
}
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

  - `damp_horizontal_wind_on_rhs::A`: Switch for applying the RHS sponge to the horizontal wind.

  - `relax_to_mean::A`: Switch for relaxing the wind towards its averages on the terrain-following surfaces. If set to `false`, the relaxation wind is computed with `relaxed_u`, `relaxed_v` and `relaxed_w`.

  - `lhs_sponge::B`: Function used to compute the Rayleigh-damping coefficient of the LHS sponge.

  - `rhs_sponge::C`: Function used to compute the Rayleigh-damping coefficient of the RHS sponge.

  - `relaxed_u::D`: Function used to compute the zonal relaxation wind if `relax_to_mean` is set to `false`.

  - `relaxed_v::E`: Function used to compute the meridional relaxation wind if `relax_to_mean` is set to `false`.

  - `relaxed_w::F`: Function used to compute the vertical relaxation wind if `relax_to_mean` is set to `false`.
"""
struct SpongeNamelist{
    A <: Bool,
    B <: Function,
    C <: Function,
    D <: Function,
    E <: Function,
    F <: Function,
}
    damp_horizontal_wind_on_rhs::A
    relax_to_mean::A
    lhs_sponge::B
    rhs_sponge::C
    relaxed_u::D
    relaxed_v::E
    relaxed_w::F
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
