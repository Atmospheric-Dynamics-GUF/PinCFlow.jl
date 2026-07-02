"""
```julia
compute_orographic_flow(
    state::State,
    deltah::AbstractFloat,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
```

Return estimates of the quantities ``\\left(\\bar{\\rho}, N^2, u_{\\mathrm{b}}, v_{\\mathrm{b}}\\right)`` of the mountain-wave-generating background flow by dispatching to the appropriate method.

```julia
compute_orographic_flow(
    orographic_flow::Val{:Surface},
    state::State,
    deltah::AbstractFloat,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
```

Return the values of ``\\left(\\bar{\\rho}, N^2, u_{\\mathrm{b}}, v_{\\mathrm{b}}\\right)`` in the first layer above the resolved surface.

```julia
compute_orographic_flow(
    orographic_flow::Val{:Summit},
    state::State,
    deltah::AbstractFloat,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
```

Return the values of ``\\left(\\bar{\\rho}, N^2, u_{\\mathrm{b}}, v_{\\mathrm{b}}\\right)`` in the first layer completely above the summit, which is assumed to be at `zctilde[i, j, k0 - 1] + deltah`.

```julia
compute_orographic_flow(
    orographic_flow::Val{:Average},
    state::State,
    deltah::AbstractFloat,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
```

Return the averages of ``\\left(\\bar{\\rho}, N^2, u_{\\mathrm{b}}, v_{\\mathrm{b}}\\right)`` between the surface and the upper edge of that layer which contains the summit.

# Arguments

  - `orographic_flow`: Method used to approximate quantities of the mountain-wave-generating background flow.

  - `state`: Model state.

  - `deltah`: Elevation difference between the local background orography and the summits of the full local orography.

  - `i`: Zonal grid index.

  - `j`: Meridional grid index.
"""
function compute_orographic_flow end

function compute_orographic_flow(
    state::State,
    deltah::AbstractFloat,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
    (; orographic_flow) = state.namelists.wkb

    @dispatch_orographic_flow return compute_orographic_flow(
        Val(orographic_flow),
        state,
        deltah,
        i,
        j,
    )
end

@ivy function compute_orographic_flow(
    orographic_flow::Val{:Surface},
    state::State,
    deltah::AbstractFloat,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
    (; k0) = state.domain
    (; rhobar, n2) = state.atmosphere
    (; u, v) = state.variables.predictands

    k = k0 - 1

    rhoh = rhobar[i, j, k]
    n2h = n2[i, j, k]
    uh = (u[i, j, k] + u[i - 1, j, k]) / 2
    vh = (v[i, j, k] + v[i, j - 1, k]) / 2

    return (rhoh, n2h, uh, vh)
end

@ivy function compute_orographic_flow(
    orographic_flow::Val{:Summit},
    state::State,
    deltah::AbstractFloat,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
    (; k0) = state.domain
    (; zctilde) = state.grid
    (; rhobar, n2) = state.atmosphere
    (; u, v) = state.variables.predictands

    k = get_next_half_level(i, j, zctilde[i, j, k0 - 1] + deltah, state) + 1

    rhoh = rhobar[i, j, k]
    n2h = n2[i, j, k]
    uh = (u[i, j, k] + u[i - 1, j, k]) / 2
    vh = (v[i, j, k] + v[i, j - 1, k]) / 2

    return (rhoh, n2h, uh, vh)
end

@ivy function compute_orographic_flow(
    orographic_flow::Val{:Average},
    state::State,
    deltah::AbstractFloat,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
    (; k0, k1) = state.domain
    (; jac, dz, zctilde) = state.grid
    (; rhobar, n2) = state.atmosphere
    (; u, v) = state.variables.predictands

    dzh = 0.0
    rhoh = 0.0
    n2h = 0.0
    uh = 0.0
    vh = 0.0
    for k in k0:k1
        dzh += jac[i, j, k] * dz
        rhoh += rhobar[i, j, k] * jac[i, j, k] * dz
        n2h += n2[i, j, k] * jac[i, j, k] * dz
        uh += (u[i, j, k] + u[i - 1, j, k]) / 2 * jac[i, j, k] * dz
        vh += (v[i, j, k] + v[i, j - 1, k]) / 2 * jac[i, j, k] * dz
        if zctilde[i, j, k] > zctilde[i, j, k0 - 1] + deltah
            break
        end
    end
    rhoh /= dzh
    n2h /= dzh
    uh /= dzh
    vh /= dzh

    return (rhoh, n2h, uh, vh)
end
