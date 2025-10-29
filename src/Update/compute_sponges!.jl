"""
```julia
compute_sponges!(state::State, dt::AbstractFloat, time::AbstractFloat)
```

Compute the Rayleigh-damping coefficients of the two sponges.

The coefficients are computed with the functions `lhs_sponge` and `rhs_sponge` in `state.namelists.sponge`.

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `time`: Simulation time.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)
"""
function compute_sponges! end

function compute_sponges!(state::State, dt::AbstractFloat, time::AbstractFloat)
    (; namelists, domain) = state
    (; lhs_sponge, rhs_sponge) = namelists.sponge
    (; zz_size, nzz, ko, i0, i1, j0, j1, k0, k1, io, jo) = domain
    (; lref, tref) = state.constants
    (; x, y, zc) = state.grid
    (; alphar, betar) = state.sponge

    kmin = ko == 0 ? k0 : k0 - 1
    kmax = ko + nzz == zz_size ? k1 : k1 + 1

    @ivy for k in kmin:kmax, j in j0:j1, i in i0:i1
        xdim = x[io + i] * lref
        ydim = y[jo + j] * lref
        zcdim = zc[i, j, k] * lref
        tdim = time * tref
        dtdim = dt * tref

        alphar[i, j, k] = lhs_sponge(xdim, ydim, zcdim, tdim, dtdim) * tref
        betar[i, j, k] = rhs_sponge(xdim, ydim, zcdim, tdim, dtdim) * tref
    end

    set_zonal_boundaries_of_field!(alphar, namelists, domain)
    set_zonal_boundaries_of_field!(betar, namelists, domain)

    set_meridional_boundaries_of_field!(alphar, namelists, domain)
    set_meridional_boundaries_of_field!(betar, namelists, domain)

    if ko == 0
        alphar[:, :, k0 - 1] .= alphar[:, :, k0]
        betar[:, :, k0 - 1] .= betar[:, :, k0]
    end
    if ko + nzz == zz_size
        alphar[:, :, k1 + 1] .= alphar[:, :, k1]
        betar[:, :, k1 + 1] .= betar[:, :, k1]
    end

    return
end
