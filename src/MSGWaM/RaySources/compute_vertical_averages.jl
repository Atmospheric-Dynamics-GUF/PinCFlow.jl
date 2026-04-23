"""
```julia
compute_vertical_averages(
    state::State,
    deltah::AbstractFloat,
    i::Integer,
    j::Integer,
)::NTuple{4, <:AbstractFloat}
```

Compute and return the local vertical averages of ``\\bar{\\rho}``, ``N^2``, ``u_{\\mathrm{b}}``, and ``v_{\\mathrm{b}}``.

# Arguments

  - `state`: Model state.

  - `deltah`: Elevation difference between the local background orography and the summits of the full local orography.

  - `i`: Zonal grid index.

  - `j`: Meridional grid index.
"""
function compute_vertical_averages end

function compute_vertical_averages(
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
    @ivy for k in k0:k1
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
