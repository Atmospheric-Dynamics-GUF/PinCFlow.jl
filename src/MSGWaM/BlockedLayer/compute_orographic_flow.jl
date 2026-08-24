"""
```julia
compute_orographic_flow(
    state::State,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
```

Return estimates of the quantities ``\\left(\\bar{\\rho}, N^2, u_{\\mathrm{b}}, v_{\\mathrm{b}}\\right)`` of the mountain-wave-generating background flow by dispatching to the appropriate method.

```julia
compute_orographic_flow(
    orographic_flow::Val{:Surface},
    state::State,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
```

Return the values of ``\\left(\\bar{\\rho}, N^2, u_{\\mathrm{b}}, v_{\\mathrm{b}}\\right)`` in the first layer above the resolved surface.

```julia
compute_orographic_flow(
    orographic_flow::Val{:Summit},
    state::State,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
```

Return the values of ``\\left(\\bar{\\rho}, N^2, u_{\\mathrm{b}}, v_{\\mathrm{b}}\\right)`` in the layer containing the summit, which is assumed to be at ``h_\\mathrm{b} + \\Delta h`` (with ``\\Delta h`` computed by `compute_elevation_difference`).

```julia
compute_orographic_flow(
    orographic_flow::Val{:Average},
    state::State,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
```

Return the averages of ``\\left(\\bar{\\rho}, N^2, u_{\\mathrm{b}}, v_{\\mathrm{b}}\\right)`` between the surface and the upper edge of that layer which contains the summit.

# Arguments

  - `orographic_flow`: Method used to approximate quantities of the mountain-wave-generating background flow.

  - `state`: Model state.

  - `i`: Zonal grid index.

  - `j`: Meridional grid index.

# See also

  - [`PinCFlow.MSGWaM.BlockedLayer.compute_elevation_difference`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.get_next_half_level`](@ref)

!!! danger "Experimental"
    The blocked-layer scheme is an experimental feature that hasn't been fully validated yet.
"""
function compute_orographic_flow end

function compute_orographic_flow(
    state::State,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
    (; orographic_flow) = state.namelists.wkb

    @dispatch_orographic_flow return compute_orographic_flow(
        Val(orographic_flow),
        state,
        i,
        j,
    )
end

@ivy function compute_orographic_flow(
    orographic_flow::Val{:Surface},
    state::State,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
    (; k0) = state.domain
    (; rhobar, n2) = state.atmosphere
    (; u, v) = state.variables.predictands

    k = k0

    rhoh = rhobar[i, j, k]
    n2h = n2[i, j, k]
    uh = (u[i, j, k] + u[i - 1, j, k]) / 2
    vh = (v[i, j, k] + v[i, j - 1, k]) / 2

    return (rhoh, n2h, uh, vh)
end

@ivy function compute_orographic_flow(
    orographic_flow::Val{:Summit},
    state::State,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
    (; k0) = state.domain
    (; hb) = state.grid
    (; rhobar, n2) = state.atmosphere
    (; u, v) = state.variables.predictands

    deltah = compute_elevation_difference(state, i, j)
    k = get_next_half_level(i, j, hb[i, j] + deltah, state)

    rhoh = rhobar[i, j, k]
    n2h = n2[i, j, k]
    uh = (u[i, j, k] + u[i - 1, j, k]) / 2
    vh = (v[i, j, k] + v[i, j - 1, k]) / 2

    return (rhoh, n2h, uh, vh)
end

@ivy function compute_orographic_flow(
    orographic_flow::Val{:Average},
    state::State,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
    (; k0, k1) = state.domain
    (; jac, dz, zctilde, hb) = state.grid
    (; rhobar, n2) = state.atmosphere
    (; u, v) = state.variables.predictands

    deltah = compute_elevation_difference(state, i, j)

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
        if zctilde[i, j, k] > hb[i, j] + deltah
            break
        end
    end
    rhoh /= dzh
    n2h /= dzh
    uh /= dzh
    vh /= dzh

    return (rhoh, n2h, uh, vh)
end
