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

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::U,
)
    (; nbx, nby, nbz) = state.namelists.domain
    (; u) = state.variables.predictands
    (; nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met, dz) = state.grid
    (; km) = state.turbulence.turbulencediffusioncoefficients

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)
    ath, bth, cth, fth, qth, pth, qth_bc, fth_bc =
        thomas_arrays(state, nx + 1, ny, nz)

    @ivy for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
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

        ith = i - nbx + 1
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

    u[(i0 - 1):i1, j0:j1, k0:k1] .= fth
    return
end

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::V,
)
    (; nbx, nby, nbz) = state.namelists.domain
    (; v) = state.variables.predictands
    (; nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met, dz) = state.grid
    (; km) = state.turbulence.turbulencediffusioncoefficients

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)
    ath, bth, cth, fth, qth, pth, qth_bc, fth_bc =
        thomas_arrays(state, nx, ny + 1, nz)

    @ivy for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
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
        jth = j - nby + 1
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

    v[i0:i1, (j0 - 1):j1, k0:k1] .= fth
    return
end

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::W,
)
    (; nbx, nby, nbz) = state.namelists.domain
    (; u, v, w) = state.variables.predictands
    (; uold, vold) = state.variables.backups
    (; nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met, dz) = state.grid
    (; km) = state.turbulence.turbulencediffusioncoefficients

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)
    ath, bth, cth, fth, qth, pth, qth_bc, fth_bc =
        thomas_arrays(state, nx, ny, nz + 1)

    @ivy for k in (k0 - 1):k1, j in j0:j1, i in i0:i1
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
        kth = k - nbz + 1

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
            2.0 * dtdz2 * met23c / jacc *
            (k33u * vu - (k33u + k33d) * vc + k33d * vd) -
            2.0 * dtdz2 / jacc^2.0 * (
                k33u * jacu * met23u * vu - (k33u + k33d) * jacc * met23c * vc +
                k33d * jacd * met23d * vd
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

    w[i0:i1, j0:j1, (k0 - 1):k1] .= fth
    return
end

function turbulent_diffusion_new!(state::State, dt::AbstractFloat, variable::W)
    (; nbx, nby, nbz) = state.namelists.domain
    (; u, v, w) = state.variables.predictands
    (; uold, vold) = state.variables.backups
    (; nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met, dz) = state.grid
    (; km) = state.turbulence.turbulencediffusioncoefficients

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)
    ath, bth, cth, fth, qth, pth, qth_bc, fth_bc =
        thomas_arrays(state, nx, ny, nz + 1)

    ii = i0:i1
    jj = j0:j1

    @ivy for k in 1:(nz + 1)
        knbz = k + nbz
        kmd =
            (
                jac[ii, jj, knbz - 1] .* km[ii, jj, knbz] .+
                jac[ii, jj, knbz] .* km[ii, jj, knbz - 1]
            ) ./ (jac[ii, jj, knbz - 1] + jac[ii, jj, knbz])
        kmu =
            (
                jac[ii, jj, knbz + 1] .* km[ii, jj, knbz] .+
                jac[ii, jj, knbz] .+ km[ii, jj, knbz + 1]
            ) ./ (jac[ii, jj, knbz + 1] .+ jac[ii, jj, knbz])
        ath[:, :, k] .= .-dtdz2 .* kmd
        bth[:, :, k] .= 1 .+ dtdz2 .* kmu .+ dtdz2 .* kmd
        cth[:, :, k] .= .-dtdz2 .* kmu

        fth[:, :, k] =
            (1 .- dtdz2 .* kmu .- dtdz2 .* kmd) .* 0.5 .*
            (w[ii, jj, knbz] .+ w[ii, jj, knbz - 1]) .+
            dtdz2 .* kmu .* 0.5 .*
            (w[ii, jj, knbz] .+ w[ii, jj, knbz + 1]) .+
            dtdz2 .* kmd .* 0.5 .* (w[ii, jj, knbz - 1] .+ w[ii, jj, knbz - 2])
    end

    thomas_algorithm!(state, ath, bth, cth, fth, qth, pth, fth_bc, qth_bc)

    w[i0:i1, j0:j1, k0:k1] .=
        0.5 .* (fth[:, :, begin:(end - 1)] .+ fth[:, :, (begin + 1):end])
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
    (; nbz) = state.namelists.domain
    (; jac, dz) = state.grid
    (; kh) = state.turbulence.turbulencediffusioncoefficients

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)
    ath, bth, cth, fth, qth, pth, qth_bc, fth_bc =
        thomas_arrays(state, nx, ny, nz)

    ii = i0:i1
    jj = j0:j1
    kk = k0:k1

    @ivy for k in 1:nz
        knbz = k + nbz
        khd =
            (
                jac[ii, jj, knbz - 1] .* kh[ii, jj, knbz] .+
                jac[ii, jj, knbz] .* kh[ii, jj, knbz - 1]
            ) ./ (jac[ii, jj, knbz - 1] + jac[ii, jj, knbz])
        khu =
            (
                jac[ii, jj, knbz + 1] .* kh[ii, jj, knbz] .+
                jac[ii, jj, knbz] .* kh[ii, jj, knbz + 1]
            ) ./ (jac[ii, jj, knbz + 1] .+ jac[ii, jj, knbz])
        ath[:, :, k] .= .-dtdz2 .* khd
        bth[:, :, k] .= 1 .+ dtdz2 .* khu .+ dtdz2 .* khd
        cth[:, :, k] .= .-dtdz2 .* khu

        fth[:, :, k] =
            (1 .- dtdz2 .* khu .- dtdz2 .* khd) .*
            (p[ii, jj, knbz] ./ (rho[ii, jj, knbz] .+ rhobar[ii, jj, knbz])) .+
            dtdz2 .* khu .* (
                p[ii, jj, knbz + 1] ./
                (rho[ii, jj, knbz + 1] .+ rhobar[ii, jj, knbz + 1])
            ) .+
            dtdz2 .* khd .* (
                p[ii, jj, knbz - 1] ./
                (rho[ii, jj, knbz - 1] .+ rhobar[ii, jj, knbz - 1])
            )
    end

    thomas_algorithm!(state, ath, bth, cth, fth, qth, pth, fth_bc, qth_bc)

    p[ii, jj, kk] .= fth .* (rho[ii, jj, kk] .+ rhobar[ii, jj, kk])
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

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)
    ath, bth, cth, fth, qth, pth, qth_bc, fth_bc =
        thomas_arrays(state, nx, ny, nz)

    ii = i0:i1
    jj = j0:j1
    kk = k0:k1

    for field in 1:fieldcount(TracerPredictands)
        chi = getfield(tracerpredictands, field)
        @ivy for k in 1:nz
            knbz = k + nbz
            khd =
                (
                    jac[ii, jj, knbz - 1] .* kh[ii, jj, knbz] .+
                    jac[ii, jj, knbz] .* kh[ii, jj, knbz - 1]
                ) ./ (jac[ii, jj, knbz - 1] + jac[ii, jj, knbz])
            khu =
                (
                    jac[ii, jj, knbz + 1] .* kh[ii, jj, knbz] .+
                    jac[ii, jj, knbz] .* kh[ii, jj, knbz + 1]
                ) ./ (jac[ii, jj, knbz + 1] .+ jac[ii, jj, knbz])
            ath[:, :, k] .= .-dtdz2 .* khd
            bth[:, :, k] .= 1 .+ dtdz2 .* khu .+ dtdz2 .* khd
            cth[:, :, k] .= .-dtdz2 .* khu

            fth[:, :, k] =
                (1 .- dtdz2 .* khu .- dtdz2 .* khd) .* chi[ii, jj, knbz] .+
                dtdz2 .* khu .* chi[ii, jj, knbz + 1] .+
                dtdz2 .* khd .* chi[ii, jj, knbz - 1]
        end

        thomas_algorithm!(state, ath, bth, cth, fth, qth, pth, fth_bc, qth_bc)

        chi[ii, jj, kk] .= fth
    end
    return
end
