function apply_saturation_scheme!(state::State, dt::AbstractFloat)
    (; testcase) = state.namelists.setting
    apply_saturation_scheme!(state, dt, testcase)
    return
end

function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    testcase::AbstractTestCase,
)
    return
end

function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    testcase::AbstractWKBTestCase,
)
    (; wkb_mode) = state.namelists.wkb
    apply_saturation_scheme!(state, dt, wkb_mode)
    return
end

function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    wkb_mode::SteadyState,
)
    return
end

function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    wkb_mode::Union{SingleColumn, MultiColumn},
)
    (; domain, grid) = state
    (; nray, rays, diffusion) = state.wkb
    (; sizex, sizey) = state.namelists.domain
    (; alpha_sat) = state.namelists.wkb
    (; io, jo, i0, i1, j0, j1, k0, k1) = state.domain
    (; lx, ly, dx, dy, ztfc) = state.grid

    for kz in k0:k1, jy in j0:j1, ix in i0:i1

        # Compute saturation integrals for wave-action reduction.
        (mb2, mb2k2) = compute_saturation_integrals(state, (ix, jy, kz))

        # Calculate the turbulent eddy diffusivity.
        n2r = interpolate_stratification(ztfc[ix, jy, kz], state, N2())
        if mb2k2 == 0 || mb2 < alpha_sat^2 * n2r^2
            diffusion[ix, jy, kz] = 0
        else
            diffusion[ix, jy, kz] =
                (mb2 - alpha_sat^2 * n2r^2) / (2 * dt * mb2k2)
        end

        # Reduce the wave-action density.
        for iray in 1:nray[ix, jy, kz]

            # Skip ray volumes with zero wave-action density.
            if rays.dens[iray, ix, jy, kz] == 0.0
                continue
            end

            xr = rays.x[iray, ix, jy, kz]
            yr = rays.y[iray, ix, jy, kz]
            zr = rays.z[iray, ix, jy, kz]

            if sizex > 1
                ix = round(Int, (xr - lx[1] - dx / 2) / dx) + i0 - io
            else
                ix = i0
            end

            if sizey > 1
                jy = round(Int, (yr - ly[1] - dy / 2) / dy) + j0 - jo
            else
                jy = j0
            end

            kz = get_next_half_level(ix, jy, zr, domain, grid)

            wnrk = rays.k[iray, ix, jy, kz]
            wnrl = rays.l[iray, ix, jy, kz]
            wnrm = rays.m[iray, ix, jy, kz]

            kappa = diffusion[ix, jy, kz]

            rays.dens[iray, ix, jy, kz] *=
                max(0, 1 - dt * 2 * kappa * (wnrk^2 + wnrl^2 + wnrm^2))
        end

        # Compute the saturation integrals again for diagnostics.
        (mb2, mb2k2) = compute_saturation_integrals(state, (ix, jy, kz))

        # Check if saturation is violated.
        n2r = interpolate_stratification(ztfc[ix, jy, kz], state, N2())
        if mb2 - alpha_sat^2 * n2r^2 > 1.0E-3 * alpha_sat^2 * n2r^2
            println("Saturation violated at (ix, jy, kz) = ", (ix, jy, kz))
            println("mb2[ix, jy, kz] = ", mb2)
            println("alpha_sat^2 * n2r^2 = ", alpha_sat^2 * n2r^2)
            println("")
        end
    end

    # Rmove rays with zero wave-action density.
    remove_rays!(state)

    return
end
