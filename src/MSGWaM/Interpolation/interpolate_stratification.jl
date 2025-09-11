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

This method first determines the two points in ``z + J \\Delta \\widehat{z} / 2`` that are closest to `zlc`. As for ``N^2``, `(i0, j0)` is used as the horizontal position, and simple linear interpolation is performed to find ``\\partial N^2 / \\partial z`` at `zlc`.

# Arguments

  - `zlc`: Vertical position of interest.

  - `state`: Model state.

  - `strtype`: Stratification quantity to interpolate.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_level`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.get_next_half_level`](@ref)
"""
function interpolate_stratification end

function interpolate_stratification(
    zlc::AbstractFloat,
    state::State,
    strtype::N2,
)::AbstractFloat
    (; domain, grid) = state
    (; bvsstrattfc) = state.atmosphere
    (; i0, j0) = domain
    (; ztfc) = grid

    kzu = get_next_level(i0, j0, zlc, domain, grid)
    kzd = kzu - 1

    @ivy zd = ztfc[i0, j0, kzd]
    @ivy zu = ztfc[i0, j0, kzu]
    @ivy strd = bvsstrattfc[i0, j0, kzd]
    @ivy stru = bvsstrattfc[i0, j0, kzu]

    if zu < zd
        error(
            "Error in interpolate_stratification (N2): zu = ",
            zu,
            " < zd = ",
            zd,
        )
    elseif zu == zd
        factor = 0.0
    elseif zlc > zu
        factor = 0.0
    elseif zlc > zd
        factor = (zu - zlc) / (zu - zd)
    else
        factor = 1.0
    end

    str = factor * strd + (1.0 - factor) * stru

    return str
end

function interpolate_stratification(
    zlc::AbstractFloat,
    state::State,
    strtype::DN2DZ,
)::AbstractFloat
    (; domain, grid) = state
    (; bvsstrattfc) = state.atmosphere
    (; i0, j0) = domain
    (; dz, ztildetfc, jac) = grid

    kzu = get_next_half_level(i0, j0, zlc, state)
    kzd = kzu - 1

    @ivy zd = ztildetfc[i0, j0, kzd]
    @ivy zu = ztildetfc[i0, j0, kzu]

    @ivy strd =
        (bvsstrattfc[i0, j0, kzd + 1] - bvsstrattfc[i0, j0, kzd]) / (
            2.0 * jac[i0, j0, kzd] * jac[i0, j0, kzd + 1] /
            (jac[i0, j0, kzd] + jac[i0, j0, kzd + 1])
        ) / dz
    @ivy stru =
        (bvsstrattfc[i0, j0, kzu + 1] - bvsstrattfc[i0, j0, kzu]) / (
            2.0 * jac[i0, j0, kzu] * jac[i0, j0, kzu + 1] /
            (jac[i0, j0, kzu] + jac[i0, j0, kzu + 1])
        ) / dz

    if zu < zd
        error(
            "Error in interpolate_stratification (DN2DZ): zu = ",
            zu,
            " < zd = ",
            zd,
        )
    elseif zu == zd
        factor = 0.0
    elseif zlc > zu
        factor = 0.0
    elseif zlc > zd
        factor = (zu - zlc) / (zu - zd)
    else
        factor = 1.0
    end

    str = factor * strd + (1.0 - factor) * stru

    return str
end
