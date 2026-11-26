"""
```julia
interpolate_stratification(
    zlc::AbstractFloat,
    state::State,
    strtype::N2,
)::AbstractFloat
```

Interpolate the squared buoyancy frequency (``N^2``) to `zlc` and return the result.

This method first determines the two points in ``z`` that are closest to `zlc`. As horizontal position, it uses `(i0, j0)`, which is arbitrary, since ``N^2`` has no horizontal dependence. Subsequently, simple linear interpolation is performed to find ``N^2`` at `zlc`.

```julia
interpolate_stratification(
    zlc::AbstractFloat,
    state::State,
    strtype::DN2DZ,
)::AbstractFloat
```

Interpolate the vertical derivative of the squared buoyancy frequency (``\\partial N^2 / \\partial z``) to `zlc` and return the result.

This method first determines the two points in ``z + J \\Delta \\hat{z} / 2`` that are closest to `zlc`. As for ``N^2``, `(i0, j0)` is used as the horizontal position, and simple linear interpolation is performed to find ``\\partial N^2 / \\partial z`` at `zlc`.

# Arguments

  - `zlc`: Vertical position of interest.

  - `state`: Model state.

  - `strtype`: Stratification quantity to interpolate.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_levels`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.get_half_levels`](@ref)
"""
function interpolate_stratification end

@ivy function interpolate_stratification(
    zlc::AbstractFloat,
    state::State,
    strtype::N2,
)::AbstractFloat
    (; domain, grid) = state
    (; n2) = state.atmosphere
    (; i0, j0) = domain
    (; zc) = grid

    (kd, ku) = get_levels(i0, j0, zlc, state)

    zd = zc[i0, j0, kd]
    zu = zc[i0, j0, ku]
    phid = n2[i0, j0, kd]
    phiu = n2[i0, j0, ku]

    interpolate((; phid, phiu, zd, zu, zlc))
end

@ivy function interpolate_stratification(
    zlc::AbstractFloat,
    state::State,
    strtype::DN2DZ,
)::AbstractFloat
    (; domain, grid) = state
    (; n2) = state.atmosphere
    (; i0, j0) = domain
    (; dz, zctilde, jac) = grid

    (kd, ku) = get_half_levels(i0, j0, zlc, state; dku = 1)

    zd = zctilde[i0, j0, kd]
    zu = zctilde[i0, j0, ku]

    phid =
        (n2[i0, j0, kd + 1] - n2[i0, j0, kd]) / (
            2.0 * jac[i0, j0, kd] * jac[i0, j0, kd + 1] /
            (jac[i0, j0, kd] + jac[i0, j0, kd + 1])
        ) / dz
    phiu =
        (n2[i0, j0, ku + 1] - n2[i0, j0, ku]) / (
            2.0 * jac[i0, j0, ku] * jac[i0, j0, ku + 1] /
            (jac[i0, j0, ku] + jac[i0, j0, ku + 1])
        ) / dz

    interpolate((; phid, phiu, zd, zu, zlc))
end
