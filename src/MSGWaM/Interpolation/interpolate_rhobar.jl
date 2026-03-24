"""
```julia
interpolate_rhobar(zlc::AbstractFloat, state::State, strtype::N2)::AbstractFloat
```

Interpolate the background density (``\\overline{\\rho}``) to `zlc` and return the result.

This method first determines the two points in ``z`` that are closest to `zlc`. As horizontal position, it uses `(i0, j0)`, which is arbitrary, since ``\\overline{\\rho}`` has no horizontal dependence. Subsequently, simple linear interpolation is performed to find ``\\overline{\\rho}`` at `zlc`.

# Arguments

  - `zlc`: Vertical position of interest.

  - `state`: Model state.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_level`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.get_next_half_level`](@ref)
"""
function interpolate_rhobar end

function interpolate_rhobar(zlc::AbstractFloat, state::State)::AbstractFloat
    (; domain, grid) = state
    (; rhobar) = state.atmosphere
    (; i0, j0) = domain
    (; zc) = grid

    ku = get_next_level(i0, j0, zlc, state; dkd = 1)
    kd = ku - 1

    @ivy zd = zc[i0, j0, kd]
    @ivy zu = zc[i0, j0, ku]
    @ivy rhod = rhobar[i0, j0, kd]
    @ivy rhou = rhobar[i0, j0, ku]

    if zu < zd
        error("Error in interpolate_rhobar (rhobar): zu = ", zu, " < zd = ", zd)
    elseif zu == zd
        factor = 0.0
    elseif zlc > zu
        factor = 0.0
    elseif zlc > zd
        factor = (zu - zlc) / (zu - zd)
    else
        factor = 1.0
    end

    rho = factor * rhod + (1.0 - factor) * rhou

    return rho
end
