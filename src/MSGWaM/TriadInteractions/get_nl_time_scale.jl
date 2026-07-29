function get_nl_time_scale end

function get_nl_time_scale(
    spec_tend::TriadTendencies,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    action_rel_tol::Float64,
)::Float64

    (; kp, m, delkp, delm) = spec_tend.spec_grid
    (; wavespectrum, col_int, action_ref) = spec_tend

    #----------------------------------------------------------
    # Minimum significant spectral-cell action:
    #
    # A_floor = action_rel_tol * A_ref
    #
    # where A_ref is the peak contained action of the weakest
    # initialized spectral mode.
    #----------------------------------------------------------

    action_floor =
        action_rel_tol * action_ref[]

    #----------------------------------------------------------
    # Maximum Boltzmann rate among dynamically active cells:
    #
    # rate = |St[kp,m]| / N[kp,m]
    #
    # A spectral cell is active when:
    #
    # N[kp,m] Δkp Δm > A_floor
    #----------------------------------------------------------

    max_rate = 0.0

    @ivy for mi in eachindex(m),
        kpi in eachindex(kp)

        was =
            wavespectrum[ii, jj, kk, kpi, mi]

        st =
            col_int[ii, jj, kk, kpi, mi]

        spectral_cell_width =
            delkp[kpi] * delm[mi]

        cell_action =
            was * spectral_cell_width

        if cell_action > action_floor &&
           !iszero(st)

            rate =
                abs(st) / was

            if rate > max_rate
                max_rate =
                    rate
            end
        end
    end

    return max_rate > 0.0 ?
        1.0 / max_rate :
        Inf
end