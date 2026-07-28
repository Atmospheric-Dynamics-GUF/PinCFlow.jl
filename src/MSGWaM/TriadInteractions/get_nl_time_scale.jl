function get_nl_time_scale end

function get_nl_time_scale(
    spec_tend::TriadTendencies,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    action_abs_tol::Float64,
    action_rel_tol::Float64,
    st_abs_tol::Float64,
)::AbstractFloat

    (; kp, m, delkp, delm) = spec_tend.spec_grid
    (; wavespectrum, col_int) = spec_tend

    @assert length(delkp) == length(kp)
    @assert length(delm) == length(m)

    #----------------------------------------------------------
    # Total wave action contained in this physical grid cell
    #
    # A_total = sum(N[kp, m] * Δkp * Δm)
    #----------------------------------------------------------

    total_action = 0.0

    @ivy for mi in eachindex(m), kpi in eachindex(kp)
        was = wavespectrum[ii, jj, kk, kpi, mi]
        spectral_cell_width = delkp[kpi] * delm[mi]

        if isfinite(was) &&
           was > 0.0 &&
           isfinite(spectral_cell_width) &&
           spectral_cell_width > 0.0

            total_action += was * spectral_cell_width
        end
    end

    if !isfinite(total_action) || total_action <= action_abs_tol
        return Inf
    end

    # A spectral cell is considered active only when its
    # contained wave action exceeds this threshold.
    action_cutoff = max(
        action_abs_tol,
        action_rel_tol * total_action,
    )

    #----------------------------------------------------------
    # Maximum Boltzmann rate among dynamically active cells
    #----------------------------------------------------------

    max_rate = 0.0

    @ivy for mi in eachindex(m), kpi in eachindex(kp)
        was = wavespectrum[ii, jj, kk, kpi, mi]
        st = col_int[ii, jj, kk, kpi, mi]
        spectral_cell_width = delkp[kpi] * delm[mi]

        if isfinite(was) &&
           isfinite(st) &&
           was > 0.0 &&
           isfinite(spectral_cell_width) &&
           spectral_cell_width > 0.0

            cell_action = was * spectral_cell_width

            if cell_action > action_cutoff &&
               abs(st) > st_abs_tol

                rate = abs(st) / was

                if isfinite(rate) && rate > max_rate
                    max_rate = rate
                end
            end
        end
    end

    return max_rate > 0.0 ? 1.0 / max_rate : Inf
end