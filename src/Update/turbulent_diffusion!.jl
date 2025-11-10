function turbulent_diffusion! end

function turbulent_diffusion!(state::State, dt::AbstractFloat)
    (; turbulence_scheme) = state.namelists.turbulence

    turbulent_diffusion!(state, dt, turbulence_scheme)
    return
end

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    turbulence_scheme::NoTurbulence,
)
    return
end

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    turbulence_scheme::TKEScheme,
)
    (; momentum_coupling, entropy_coupling) = state.namelists.turbulence
    if momentum_coupling
        turbulent_diffusion!(state, dt, U())
        turbulent_diffusion!(state, dt, V())
        turbulent_diffusion!(state, dt, W())
    end
    if entropy_coupling
        turbulent_diffusion!(state, dt, Theta())
    end
    turbulent_diffusion!(state, dt, Chi())
    return
end

function turbulent_diffusion!(state::State, dt::AbstractFloat, variable::U)
    (; nbx, nby, nbz) = state.namelists.domain
    (; u) = state.variables.predictands
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met, dz) = state.grid
    (; ath, bth, cth, fth, qth, pth, qth_bc, fth_bc) =
        state.variables.auxiliaries
    (; km) = state.turbulence.turbulencediffusioncoefficients
    dtdz2 = dt / (2.0 * dz^2.0)
    reset_thomas!(state)
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        k33u =
            0.5 * (
                (
                    jac[i, j, k] * (
                        jac[i, j, k + 1] *
                        met[i, j, k + 1, 3, 3] *
                        km[i, j, k + 1]
                    ) +
                    jac[i, j, k + 1] *
                    (jac[i, j, k] * met[i, j, k, 3, 3] * km[i, j, k])
                ) / (jac[i, j, k] + jac[i, j, k + 1]) +
                (
                    jac[i + 1, j, k] * (
                        jac[i + 1, j, k + 1] *
                        met[i + 1, j, k + 1, 3, 3] *
                        km[i + 1, j, k + 1]
                    ) +
                    jac[i + 1, j, k + 1] * (
                        jac[i + 1, j, k] *
                        met[i + 1, j, k, 3, 3] *
                        km[i + 1, j, k]
                    )
                ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
            )

        k33d =
            0.5 * (
                (
                    jac[i, j, k] * (
                        jac[i, j, k - 1] *
                        met[i, j, k - 1, 3, 3] *
                        km[i, j, k - 1]
                    ) +
                    jac[i, j, k - 1] *
                    (jac[i, j, k] * met[i, j, k, 3, 3] * km[i, j, k])
                ) / (jac[i, j, k] + jac[i, j, k - 1]) +
                (
                    jac[i + 1, j, k] * (
                        jac[i + 1, j, k - 1] *
                        met[i + 1, j, k - 1, 3, 3] *
                        km[i + 1, j, k - 1]
                    ) +
                    jac[i + 1, j, k - 1] * (
                        jac[i + 1, j, k] *
                        met[i + 1, j, k, 3, 3] *
                        km[i + 1, j, k]
                    )
                ) / (jac[i + 1, j, k] + jac[i + 1, j, k - 1])
            )

        jacc = 0.5 * (jac[i, j, k] + jac[i + 1, j, k])

        ith = i - nbx
        jth = j - nby
        kth = k - nbz

        ath[ith, jth, kth] = -dtdz2 * k33d / jacc
        bth[ith, jth, kth] = 1 + dtdz2 * k33u / jacc + dtdz2 * k33d / jacc
        cth[ith, jth, kth] = -dtdz2 * k33u / jacc
        fth[ith, jth, kth] =
            dtdz2 * k33u / jacc * u[i, j, k + 1] +
            (1 - dtdz2 * k33u / jacc - dtdz2 * k33d / jacc) * u[i, j, k] +
            dtdz2 * k33d / jacc * u[i, j, k - 1]
    end

    thomas_algorithm!(state, ath, bth, cth, fth, qth, pth, fth_bc, qth_bc)

    u[i0:i1, j0:j1, k0:k1] .= fth
    return
end

function turbulent_diffusion!(state::State, dt::AbstractFloat, variable::V)
    (; nbx, nby, nbz) = state.namelists.domain
    (; v) = state.variables.predictands
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met, dz) = state.grid
    (; ath, bth, cth, fth, qth, pth, qth_bc, fth_bc) =
        state.variables.auxiliaries
    (; km) = state.turbulence.turbulencediffusioncoefficients
    dtdz2 = dt / (2.0 * dz^2.0)
    reset_thomas!(state)
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        k33u =
            0.5 * (
                (
                    jac[i, j, k] * (
                        jac[i, j, k + 1] *
                        met[i, j, k + 1, 3, 3] *
                        km[i, j, k + 1]
                    ) +
                    jac[i, j, k + 1] *
                    (jac[i, j, k] * met[i, j, k, 3, 3] * km[i, j, k])
                ) / (jac[i, j, k] + jac[i, j, k + 1]) +
                (
                    jac[i, j + 1, k] * (
                        jac[i, j + 1, k + 1] *
                        met[i, j + 1, k + 1, 3, 3] *
                        km[i, j + 1, k + 1]
                    ) +
                    jac[i, j + 1, k + 1] * (
                        jac[i, j + 1, k] *
                        met[i, j + 1, k, 3, 3] *
                        km[i, j + 1, k]
                    )
                ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
            )

        k33d =
            0.5 * (
                (
                    jac[i, j, k] * (
                        jac[i, j, k - 1] *
                        met[i, j, k - 1, 3, 3] *
                        km[i, j, k - 1]
                    ) +
                    jac[i, j, k - 1] *
                    (jac[i, j, k] * met[i, j, k, 3, 3] * km[i, j, k])
                ) / (jac[i, j, k] + jac[i, j, k - 1]) +
                (
                    jac[i, j + 1, k] * (
                        jac[i, j + 1, k - 1] *
                        met[i, j + 1, k - 1, 3, 3] *
                        km[i, j + 1, k - 1]
                    ) +
                    jac[i, j + 1, k - 1] * (
                        jac[i, j + 1, k] *
                        met[i, j + 1, k, 3, 3] *
                        km[i, j + 1, k]
                    )
                ) / (jac[i, j + 1, k] + jac[i, j + 1, k - 1])
            )

        jacc = 0.5 * (jac[i, j, k] + jac[i, j + 1, k])

        ith = i - nbx
        jth = j - nby
        kth = k - nbz

        ath[ith, jth, kth] = -dtdz2 * k33d / jacc
        bth[ith, jth, kth] = 1 + dtdz2 * k33u / jacc + dtdz2 * k33d / jacc
        cth[ith, jth, kth] = -dtdz2 * k33u / jacc
        fth[ith, jth, kth] =
            dtdz2 * k33u / jacc * v[i, j, k + 1] +
            (1 - dtdz2 * k33u / jacc - dtdz2 * k33d / jacc) * v[i, j, k] +
            dtdz2 * k33d / jacc * v[i, j, k - 1]
    end

    thomas_algorithm!(state, ath, bth, cth, fth, qth, pth, fth_bc, qth_bc)

    v[i0:i1, j0:j1, k0:k1] .= fth
    return
end

function turbulent_diffusion!(state::State, dt::AbstractFloat, variable::W)
    (; nbx, nby, nbz) = state.namelists.domain
    (; u, v, w) = state.variables.predictands
    (; uold, vold) = state.variables.backups
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met, dz) = state.grid
    (; ath, bth, cth, fth, qth, pth, qth_bc, fth_bc) =
        state.variables.auxiliaries
    (; km) = state.turbulence.turbulencediffusioncoefficients
    dtdz2 = dt / (2.0 * dz^2.0)
    reset_thomas!(state)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        met13c =
            (
                jac[i, j, k] * met[i, j, k + 1, 1, 3] +
                jac[i, j, k + 1] * met[i, j, k, 1, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met23c =
            (
                jac[i, j, k] * met[i, j, k + 1, 2, 3] +
                jac[i, j, k + 1] * met[i, j, k, 2, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        jacc =
            2.0 * jac[i, j, k] * jac[i, j, k + 1] /
            (jac[i, j, k] + jac[i, j, k + 1])
        met13u =
            (
                jac[i, j, k + 1] * met[i, j, k + 2, 1, 3] +
                jac[i, j, k + 2] * met[i, j, k + 1, 1, 3]
            ) / (jac[i, j, k + 1] + jac[i, j, k + 2])
        met23u =
            (
                jac[i, j, k + 1] * met[i, j, k + 2, 2, 3] +
                jac[i, j, k + 2] * met[i, j, k + 1, 2, 3]
            ) / (jac[i, j, k + 1] + jac[i, j, k + 2])
        jacu =
            2.0 * jac[i, j, k + 1] * jac[i, j, k + 2] /
            (jac[i, j, k + 1] + jac[i, j, k + 2])
        met13d =
            (
                jac[i, j, k] * met[i, j, k - 1, 1, 3] +
                jac[i, j, k - 1] * met[i, j, k, 1, 3]
            ) / (jac[i, j, k] + jac[i, j, k - 1])
        met23d =
            (
                jac[i, j, k] * met[i, j, k - 1, 2, 3] +
                jac[i, j, k - 1] * met[i, j, k, 2, 3]
            ) / (jac[i, j, k] + jac[i, j, k - 1])
        jacd =
            2.0 * jac[i, j, k] * jac[i, j, k - 1] /
            (jac[i, j, k] + jac[i, j, k - 1])
        k33d = jac[i, j, k] * met[i, j, k, 3, 3] * km[i, j, k]
        k33u = jac[i, j, k + 1] * met[i, j, k + 1, 3, 3] * km[i, j, k + 1]

        uc =
            (
                jac[i, j, k] *
                0.25 *
                (
                    u[i - 1, j, k + 1] +
                    u[i, j, k + 1] +
                    uold[i - 1, j, k + 1] +
                    uold[i, j, k + 1]
                ) +
                jac[i, j, k + 1] *
                0.25 *
                (
                    u[i - 1, j, k] +
                    u[i, j, k] +
                    uold[i - 1, j, k] +
                    uold[i, j, k]
                )
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        vc =
            (
                jac[i, j, k] *
                0.25 *
                (
                    v[i, j - 1, k + 1] +
                    v[i, j, k + 1] +
                    vold[i, j - 1, k + 1] +
                    vold[i, j, k + 1]
                ) +
                jac[i, j, k + 1] *
                0.25 *
                (
                    v[i, j - 1, k] +
                    v[i, j, k] +
                    vold[i, j - 1, k] +
                    vold[i, j, k]
                )
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        ud =
            (
                jac[i, j, k - 1] *
                0.25 *
                (
                    u[i - 1, j, k] +
                    u[i, j, k] +
                    uold[i - 1, j, k] +
                    uold[i, j, k]
                ) +
                jac[i, j, k] *
                0.25 *
                (
                    u[i - 1, j, k - 1] +
                    u[i, j, k - 1] +
                    uold[i - 1, j, k - 1] +
                    uold[i, j, k - 1]
                )
            ) / (jac[i, j, k - 1] + jac[i, j, k])
        vd =
            (
                jac[i, j, k - 1] *
                0.25 *
                (
                    v[i, j - 1, k + 1] +
                    v[i, j, k] +
                    vold[i, j - 1, k + 1] +
                    vold[i, j, k]
                ) +
                jac[i, j, k] *
                0.25 *
                (
                    v[i, j - 1, k - 1] +
                    v[i, j, k - 1] +
                    vold[i, j - 1, k - 1] +
                    vold[i, j, k - 1]
                )
            ) / (jac[i, j, k - 1] + jac[i, j, k])
        uu =
            (
                jac[i, j, k + 1] *
                0.25 *
                (
                    u[i - 1, j, k + 2] +
                    u[i, j, k + 2] +
                    uold[i - 1, j, k + 2] +
                    uold[i, j, k + 2]
                ) +
                jac[i, j, k + 2] *
                0.25 *
                (
                    u[i - 1, j, k + 1] +
                    u[i, j, k + 1] +
                    uold[i - 1, j, k + 1] +
                    uold[i, j, k + 1]
                )
            ) / (jac[i, j, k + 1] + jac[i, j, k + 2])
        vu =
            (
                jac[i, j, k + 1] *
                0.25 *
                (
                    v[i, j - 1, k + 2] +
                    v[i, j, k + 2] +
                    vold[i, j - 1, k + 2] +
                    vold[i, j, k + 2]
                ) +
                jac[i, j, k + 2] *
                0.25 *
                (
                    v[i, j - 1, k + 1] +
                    v[i, j, k + 1] +
                    vold[i, j - 1, k + 1] +
                    vold[i, j, k + 1]
                )
            ) / (jac[i, j, k + 1] + jac[i, j, k + 2])

        ith = i - nbx
        jth = j - nby
        kth = k - nbz

        ath[ith, jth, kth] = -dtdz2 / jacc^2.0 * k33d * jacd
        bth[ith, jth, kth] = 1.0 + dtdz2 / jacc^2.0 * (k33u + k33d) * jacc
        cth[ith, jth, kth] = -dtdz2 / jacc^2.0 * k33u * jacu
        fth[ith, jth, kth] =
            2.0 * dtdz2 * met13c / jacc *
            (k33u * uu - (k33u + k33d) * uc + k33d * ud) -
            2.0 * dtdz2 / jacc^2.0 * (
                k33u * jacu * met13u * uu - (k33u + k33d) * jacc * met13c * uc +
                k33d * jacd * met13d * ud
            ) +
            2.0 * dtdz2 * met13c / jacc *
            (k33u * vu - (k33u + k33d) * vc + k33d * vd) -
            2.0 * dtdz2 / jacc^2.0 * (
                k33u * jacu * met13u * vu - (k33u + k33d) * jacc * met13c * vc +
                k33d * jacd * met13d * vd
            ) +
            dtdz2 *
            jacc^2.0 *
            (
                k33u * jacu * w[i, j, k + 1] -
                (k33u + k33d) * jacc * w[i, j, k] +
                k33d * jacd * w[i, j, k - 1]
            )
    end

    thomas_algorithm!(state, ath, bth, cth, fth, qth, pth, fth_bc, qth_bc)

    w[i0:i1, j0:j1, k0:k1] .= fth
    return
end

function turbulent_diffusion!(state::State, dt::AbstractFloat, variable::Theta)
    (; model) = state.namelists.atmosphere

    turbulent_diffusion!(state, dt, variable, model)
    return
end

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::Theta,
    model::Union{PseudoIncompressible, Boussinesq},
)
    return
end

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::Theta,
    model::Compressible,
)
    return
end

function turbulent_diffusion!(state::State, dt::AbstractFloat, variable::Chi)
    (; tracer_setup) = state.namelists.tracer

    turbulent_diffusion!(state, dt, variable, tracer_setup)
    return
end

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::Chi,
    tracer_setup::NoTracer,
)
    return
end

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::Theta,
    tracer_setup::TracerOn,
)
    return
end