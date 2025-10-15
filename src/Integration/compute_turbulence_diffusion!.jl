function compute_turbulence_diffusion!(state::State)
    (; turbulence_scheme) = state.namelists.turbulence

    compute_turbulence_diffusion!(state, turbulence_scheme)
    return
end

function compute_turbulence_diffusion!(
    state::State,
    turbulence_scheme::NoTurbulence,
)
    return
end

function compute_turbulence_diffusion!(
    state::State,
    turbulence_scheme::TKEScheme,
)
    (; lturb, ck, c3) = state.turbulence.turbulenceconstants
    (; kh, km, kek) = state.turbulence.turbulenceauxiliaries
    (; tke) = state.turbulence.turbulencepredictands
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere
    (; k0, k1, j0, j1, i0, i1) = state.domain
    (; lref, tref) = state.constants

    lturb_ndim = lturb / lref

    check_tke!(state)
    kh[i0:i1, j0:j1, k0:k1] .=
        c3 .* ck .* lturb_ndim .*
        sqrt.(
            tke[i0:i1, j0:j1, k0:k1] ./
            (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )
    km[i0:i1, j0:j1, k0:k1] .=
        ck .* lturb_ndim .*
        sqrt.(
            tke[i0:i1, j0:j1, k0:k1] ./
            (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )
    kek[i0:i1, j0:j1, k0:k1] .=
        lturb_ndim .*
        sqrt.(
            tke[i0:i1, j0:j1, k0:k1] ./
            (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )

    set_boundaries!(state, BoundaryDiffusionCoefficients())
    return
end

function compute_turbulence_diffusion!(
    state::State,
    turbulence_scheme::TTEScheme,
)
    (; lturb) = state.turbulence.turbulenceconstants
    (; kh, km, kek) = state.turbulence.turbulenceauxiliaries
    (; tke) = state.turbulence.turbulencepredictands
    (; rho) = state.variables.predictands
    (; k0, k1, j0, j1, i0, i1) = state.domain

    lturb_ndim = lturb / lref

    kh[i0:i1, j0:j1, k0:k1] .=
        lturb_ndim .*
        sqrt.(
            tke[i0:i1, j0:j1, k0:k1] ./
            (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )
    km[i0:i1, j0:j1, k0:k1] .=
        lturb_ndim .*
        sqrt.(
            tke[i0:i1, j0:j1, k0:k1] ./
            (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )
    kek[i0:i1, j0:j1, k0:k1] .=
        lturb_ndim .*
        sqrt.(
            tke[i0:i1, j0:j1, k0:k1] ./
            (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )

    set_boundaries!(state, BoundaryDiffusionCoefficients())
    return
end