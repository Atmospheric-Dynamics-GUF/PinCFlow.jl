function get_dephasing_time end

function get_dephasing_time(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    tau_nl::Float64,
    ::Triad2D,
)::Float64

    (; x_size) = state.namelists.domain
    (; action_rel_tol, increment_rel_tol) = state.namelists.triad

    (; spec_tend) = state
    (; wavespectrum, col_int, diag_dephasing_time, action_ref) = spec_tend
    (; delkp, delm) = spec_tend.spec_grid

    # No action-active spectral cell has a nonzero collision integral.
    if isinf(tau_nl)
        return Inf
    end

    action_floor = action_rel_tol * action_ref[]

    # tau_nl is defined from the maximum active |St| / N.
    max_rate = 1.0 / tau_nl
    rate_cutoff = increment_rel_tol * max_rate

    dephasing_time_min = Inf

    @ivy for mi in axes(diag_dephasing_time, 2), kpi in axes(diag_dephasing_time, 1)

        tau_dep = diag_dephasing_time[kpi, mi]

        # No contributing resonant triad associated with this target cell,
        # or the relevant triads do not dephase.
        if !isfinite(tau_dep)
            continue
        end

        was = wavespectrum[ii, jj, kk, kpi, mi]

        spectral_cell_measure = x_size == 1 ? abs(delm[mi]) : abs(delkp[kpi] * delm[mi])
        cell_action = was * spectral_cell_measure

        # Ignore target spectral cells with negligible contained action.
        if cell_action <= action_floor
            continue
        end

        # was > 0 because cell_action > action_floor.
        st = col_int[ii, jj, kk, kpi, mi]
        rate = abs(st) / was

        # Ignore target cells whose net nonlinear tendency is negligible.
        if rate <= rate_cutoff
            continue
        end

        dephasing_time_min = min(dephasing_time_min, tau_dep)
    end

    return dephasing_time_min
end