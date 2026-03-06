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
    (; momentum_coupling, entropy_coupling, tracer_coupling) =
        state.namelists.turbulence
    (; uold, vold) = state.variables.backups
    (; u, v) = state.variables.predictands

    if momentum_coupling
        uold .= copy(u)
        vold .= copy(v)
        turbulent_diffusion!(state, dt, U())
        turbulent_diffusion!(state, dt, V())
        turbulent_diffusion!(state, dt, W())
    end
    if entropy_coupling
        turbulent_diffusion!(state, dt, Theta())
    end
    if tracer_coupling
        turbulent_diffusion!(state, dt, Chi())
    end
    return
end

function turbulent_diffusion!(state::State, dt::AbstractFloat, variable::U)
    (; nbx, nby, nbz) = state.namelists.domain
    (; u) = state.variables.predictands
    (; nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met, dz) = state.grid
    (; km) = state.turbulence.turbulencediffusioncoefficients
    (; ath, bth, cth, fth, qth, pth, qth_bc, fth_bc) =
        state.variables.auxiliaries

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
    (; nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met, dz) = state.grid
    (; km) = state.turbulence.turbulencediffusioncoefficients
    (; ath, bth, cth, fth, qth, pth, qth_bc, fth_bc) =
        state.variables.auxiliaries

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
    (; nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met, dz) = state.grid
    (; km) = state.turbulence.turbulencediffusioncoefficients
    (; ath, bth, cth, fth, qth, pth, qth_bc, fth_bc) =
        state.variables.auxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        k33u = jac[i, j, k + 1] * met[i, j, k + 1, 3, 3] * km[i, j, k + 1]
        k33d = jac[i, j, k - 1] * met[i, j, k - 1, 3, 3] * km[i, j, k - 1]

        wu = compute_vertical_wind(i, j, k + 1, state)
        wc = compute_vertical_wind(i, j, k, state)
        wd = compute_vertical_wind(i, j, k - 1, state)

        dudtd =
            0.5 * (
                (u[i - 1, j, k] - uold[i - 1, j, k]) / dt +
                (u[i, j, k] - uold[i, j, k]) / dt
            )
        dudtu =
            0.5 * (
                (u[i - 1, j, k + 1] - uold[i - 1, j, k + 1]) / dt +
                (u[i, j, k + 1] - uold[i, j, k + 1]) / dt
            )
        dudtc13 =
            (
                jac[i, j, k] * met[i, j, k + 1, 1, 3] * dudtu +
                jac[i, j, k + 1] * met[i, j, k, 1, 3] * dudtd
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        dvdtd =
            0.5 * (
                (v[i, j - 1, k] - vold[i, j - 1, k]) / dt +
                (v[i, j, k] - vold[i, j, k]) / dt
            )
        dvdtu =
            0.5 * (
                (v[i, j - 1, k + 1] - vold[i, j - 1, k + 1]) / dt +
                (v[i, j, k + 1] - vold[i, j, k + 1]) / dt
            )
        dvdtc23 =
            (
                jac[i, j, k] * met[i, j, k + 1, 2, 3] * dvdtu +
                jac[i, j, k + 1] * met[i, j, k, 2, 3] * dvdtd
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        jacc =
            2.0 * jac[i, j, k] * jac[i, j, k + 1] /
            (jac[i, j, k] + jac[i, j, k + 1])

        ith = i - nbx
        jth = j - nby
        kth = k - nbz

        ath[ith, jth, kth] = -dtdz2 / jacc * k33d
        bth[ith, jth, kth] = 1 + dtdz2 / jacc * k33u + dtdz2 / jacc * k33d
        cth[ith, jth, kth] = -dtdz2 / jacc * k33u
        fth[ith, jth, kth] =
            dtdz2 / jacc * k33u * wu +
            (1 - dtdz2 / jacc * k33u - dtdz2 / jacc * k33d) * wc +
            dtdz2 / jacc * k33d * wd +
            dudtc13 +
            dvdtc23
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
    (; p, rho) = state.variables.predictands
    (; rhobar) = state.atmosphere
    (; nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (; nbx, nby, nbz) = state.namelists.domain
    (; jac, dz) = state.grid
    (; kh) = state.turbulence.turbulencediffusioncoefficients
    (; ath, bth, cth, fth, qth, pth, qth_bc, fth_bc) =
        state.variables.auxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        khd =
            (jac[i, j, k - 1] * kh[i, j, k] + jac[i, j, k] * kh[i, j, k - 1]) /
            (jac[i, j, k - 1] + jac[i, j, k])
        khu =
            (jac[i, j, k + 1] * kh[i, j, k] + jac[i, j, k] * kh[i, j, k + 1]) /
            (jac[i, j, k + 1] + jac[i, j, k])

        ith = i - nbx
        jth = j - nby
        kth = k - nbz

        ath[ith, jth, kth] = -dtdz2 * khd
        bth[ith, jth, kth] = 1 + dtdz2 * khu + dtdz2 * khd
        cth[ith, jth, kth] = -dtdz2 * khu

        fth[ith, jth, kth] =
            (1 - dtdz2 * khu - dtdz2 * khd) * p[i, j, k] +
            #(p[i, j, k] / (rho[i, j, k] + rhobar[i, j, k])) +
            dtdz2 * khu * p[i, j, k + 1] +
            #(p[i, j, k + 1] / (rho[i, j, k + 1] + rhobar[i, j, k + 1])) +
            dtdz2 * khd * p[i, j, k - 1]
        #(p[i, j, k - 1] / (rho[i, j, k - 1] + rhobar[i, j, k - 1]))
    end

    thomas_algorithm!(state, ath, bth, cth, fth, qth, pth, fth_bc, qth_bc)

    p[i0:i1, j0:j1, k0:k1] .= fth #.* (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
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
    variable::Chi,
    tracer_setup::TracerOn,
)
    (; tracerpredictands) = state.tracer
    (; nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (; nbz) = state.namelists.domain
    (; jac, dz) = state.grid
    (; kh) = state.turbulence.turbulencediffusioncoefficients
    (; ath, bth, cth, fth, qth, pth, qth_bc, fth_bc) =
        state.variables.auxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    ii = i0:i1
    jj = j0:j1
    kk = k0:k1

    for field in 1:fieldcount(TracerPredictands)
        chi = getfield(tracerpredictands, field)
        @ivy for k in 1:nz
            knbz = k + nbz
            khd =
                (
                    jac[i, j, k - 1] .* kh[i, j, k] .+
                    jac[i, j, k] .* kh[i, j, k - 1]
                ) ./ (jac[i, j, k - 1] + jac[i, j, k])
            khu =
                (
                    jac[i, j, k + 1] .* kh[i, j, k] .+
                    jac[i, j, k] .* kh[i, j, k + 1]
                ) ./ (jac[i, j, k + 1] .+ jac[i, j, k])
            ath[:, :, k] .= .-dtdz2 .* khd
            bth[:, :, k] .= 1 .+ dtdz2 .* khu .+ dtdz2 .* khd
            cth[:, :, k] .= .-dtdz2 .* khu

            fth[:, :, k] =
                (1 .- dtdz2 .* khu .- dtdz2 .* khd) .* chi[i, j, k] .+
                dtdz2 .* khu .* chi[i, j, k + 1] .+
                dtdz2 .* khd .* chi[i, j, k - 1]
        end

        thomas_algorithm!(state, ath, bth, cth, fth, qth, pth, fth_bc, qth_bc)

        chi[ii, jj, kk] .= fth
    end
    return
end
