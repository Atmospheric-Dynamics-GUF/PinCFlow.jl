function compute_tke_terms!(state::State, p0::Predictands)
    (; dissipation, diffusion, advection) =
        state.turbulence.turbulenceauxiliaries
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; tke) = state.turbulence.turbulencepredictands
    (; lturb_ndim) = state.turbulence.turbulenceconstants
    (; rhobar) = state.atmosphere
    (; dx, dy, dz, jac) = state.grid
    (; rho) = p0
    (; kek) = state.turbulence.turbulencediffusioncoefficients
    (; phitke) = state.turbulence.turbulencefluxes
    (; turbulence_scheme) = state.namelists.turbulence

    reconstruct!(state, turbulence_scheme)
    set_boundaries!(state, BoundaryReconstructions())
    compute_fluxes!(state, p0, turbulence_scheme)
    set_boundaries!(state, BoundaryFluxes())

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        specific_tke = tke[i, j, k] / (rho[i, j, k] + rhobar[i, j, k])
        dissipation[i, j, k] =
            -2 * tke[i, j, k] * sqrt(specific_tke) / lturb_ndim

        fl = phitke[i - 1, j, k, 1]
        fr = phitke[i, j, k, 1]
        gb = phitke[i, j - 1, k, 2]
        gf = phitke[i, j, k, 2]
        hd = phitke[i, j, k - 1, 3]
        hu = phitke[i, j, k, 3]

        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
        fluxdiff /= jac[i, j, k]

        advection[i, j, k] = -fluxdiff

        kekd =
            (
                jac[i, j, k - 1] * kek[i, j, k] +
                jac[i, j, k] * kek[i, j, k - 1]
            ) / (jac[i, j, k - 1] + jac[i, j, k])
        keku =
            (
                jac[i, j, k + 1] * kek[i, j, k] +
                jac[i, j, k] * kek[i, j, k - 1]
            ) / (jac[i, j, k + 1] + jac[i, j, k])

        dekdzu = (tke[i, j, k + 1] - tke[i, j, k]) / dz 
        dekdzd = (tke[i, j, k] - tke[i, j, k - 1]) / dz

        diffusion[i, j, k] = (keku * dekdzu - kekd * dekdzd) / dz
    end

    return
end