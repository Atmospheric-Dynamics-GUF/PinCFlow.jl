function apply_spectral_dissipation! end

function apply_spectral_dissipation!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    dtau::AbstractFloat,
    ::Triad2D,
)
    (; spec_tend) = state
    (; wavespectrum) = spec_tend
    (; kp, m) = spec_tend.spec_grid
    (; n2) = state.atmosphere
    (; dissipation_strength) = state.namelists.triad

    # ----------------------------------------------------------
    # Dissipative scales are determined directly from the
    # resolved spectral grid.
    # ----------------------------------------------------------

    kd_inf = kp[1]
    kd_sup = kp[end]

    md_inf = minimum(abs, m)
    md_sup = maximum(abs, m)

    n_local = sqrt(n2[ii, jj, kk])

    @ivy for mi in eachindex(m), kpi in eachindex(kp)

        kpi_val = kp[kpi]
        mi_abs = abs(m[mi])

        if mi_abs == 0.0
            continue
        end

        # Hydrostatic intrinsic frequency.
        omega_hat = n_local * kpi_val / mi_abs

        dissipation_shape =
            (kd_inf / kpi_val)^8 +
            (md_inf / mi_abs)^8 +
            (kpi_val / kd_sup)^8 +
            (mi_abs / md_sup)^8

        dissipation_rate =
            dissipation_strength * dissipation_shape / omega_hat

        # Implicit dissipative update:
        #
        #     N^{n+1} = N^n / (1 + D Δτ)
        #
        # dtau supplied to this function may be a half timestep.
        wavespectrum[ii, jj, kk, kpi, mi] /=
            1.0 + dtau * dissipation_rate
    end

    return nothing
end