function get_dephasing_time end

function get_dephasing_time(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    tau_nl::Float64,
    ::Triad2D,
)::Float64

    (; spec_tend) = state.wkb
    (; wavespectrum, col_int, action_ref) = spec_tend
    (; kp, m, kpc, delkp, delm) = spec_tend.spec_grid
    (; aa, la, qq, lq) = spec_tend.kin_box
    (; action_rel_tol, increment_rel_tol) = state.namelists.triad

    (; n2) = state.atmosphere
    (; u) = state.variables.predictands

    # No action-active spectral cell has a nonzero collision
    # integral in this physical grid cell.
    if isinf(tau_nl)
        return Inf
    end

    # Minimum significant spectral-cell action:
    #
    # A_floor = action_rel_tol * A_ref
    action_floor = action_rel_tol * action_ref[]

    # Maximum nonlinear rate among action-active spectral cells:
    #
    # r_max = max(|St| / N) = 1 / tau_nl
    max_rate = 1.0 / tau_nl

    # Retain a parent spectral cell when:
    #
    # (|St| / N) / r_max > increment_rel_tol
    rate_cutoff = increment_rel_tol * max_rate

    n_local = sqrt(n2[ii, jj, kk])
    dudz = compute_dphidz_center(u, state, ii, jj, kk, identity)
    dndz = compute_dphidz_center(n2, state, ii, jj, kk, sqrt)

    dephasing_time_min = Inf
    eps_denom = 1.0e-14

    @ivy for mi in eachindex(m), kpi in eachindex(kp)
        was = wavespectrum[ii, jj, kk, kpi, mi]
        st = col_int[ii, jj, kk, kpi, mi]

        spectral_cell_width = delkp[kpi] * delm[mi]
        cell_action = was * spectral_cell_width

        # Exclude spectral cells whose contained action is below
        # the initialized action floor.
        if cell_action <= action_floor
            continue
        end

        rate = abs(st) / was

        # Exclude parent cells whose nonlinear rate is negligible
        # relative to the maximum active nonlinear rate.
        if rate <= rate_cutoff
            continue
        end

        kp_parent = kp[kpi]
        m_parent = m[mi]

        if abs(m_parent) < eps_denom
            continue
        end

        aar = aa[kpi]
        qqr = qq[kpi]

        #------------------------------------------------------
        # Sum interaction manifold
        #------------------------------------------------------
        if kp_parent > 2.0 * kpc[1]
            for ai in 1:la[kpi]

                #--------------------------------------------------
                # Left branch
                #--------------------------------------------------

                p_left = aar[ai] - kp_parent
                kp_1, kp_2 = compute_kp1kp2(kp_parent, p_left, Sum())

                m_1, m_2 = compute_m1m2(kp_parent, kp_1, kp_2, m_parent, Sum(), Sum())

                dephasing_time_min = update_dephasing_time(
                    dephasing_time_min,
                    n_local,
                    dudz,
                    dndz,
                    kp_parent,
                    m_parent,
                    kp_1,
                    m_1,
                    kp_2,
                    m_2,
                    eps_denom,
                    Sum(),
                )

                m_1, m_2 = compute_m1m2(kp_parent, kp_1, kp_2, m_parent, Sum(), Difference())

                dephasing_time_min = update_dephasing_time(
                    dephasing_time_min,
                    n_local,
                    dudz,
                    dndz,
                    kp_parent,
                    m_parent,
                    kp_1,
                    m_1,
                    kp_2,
                    m_2,
                    eps_denom,
                    Sum(),
                )

                #--------------------------------------------------
                # Right branch
                #--------------------------------------------------

                p_right = kp_parent - aar[ai]
                kp_1, kp_2 = compute_kp1kp2(kp_parent, p_right, Sum())

                m_1, m_2 = compute_m1m2(kp_parent, kp_1, kp_2, m_parent, Sum(), Sum())

                dephasing_time_min = update_dephasing_time(
                    dephasing_time_min,
                    n_local,
                    dudz,
                    dndz,
                    kp_parent,
                    m_parent,
                    kp_1,
                    m_1,
                    kp_2,
                    m_2,
                    eps_denom,
                    Sum(),
                )

                m_1, m_2 = compute_m1m2(kp_parent, kp_1, kp_2, m_parent, Sum(), Difference())

                dephasing_time_min = update_dephasing_time(
                    dephasing_time_min,
                    n_local,
                    dudz,
                    dndz,
                    kp_parent,
                    m_parent,
                    kp_1,
                    m_1,
                    kp_2,
                    m_2,
                    eps_denom,
                    Sum(),
                )
            end
        end
        #------------------------------------------------------
        # Difference interaction manifold
        #------------------------------------------------------
        if kp_parent < kpc[end] - kpc[1]
            for qi in 1:lq[kpi]
                q_val = qqr[qi]
                kp_1, kp_2 = compute_kp1kp2(kp_parent, q_val, Difference())

                m_1, m_2 = compute_m1m2(kp_parent, kp_1, kp_2, m_parent, Difference(), Sum())

                dephasing_time_min = update_dephasing_time(
                    dephasing_time_min,
                    n_local,
                    dudz,
                    dndz,
                    kp_parent,
                    m_parent,
                    kp_1,
                    m_1,
                    kp_2,
                    m_2,
                    eps_denom,
                    Difference(),
                )

                m_1, m_2 = compute_m1m2(kp_parent, kp_1, kp_2, m_parent, Difference(), Difference())

                dephasing_time_min = update_dephasing_time(
                    dephasing_time_min,
                    n_local,
                    dudz,
                    dndz,
                    kp_parent,
                    m_parent,
                    kp_1,
                    m_1,
                    kp_2,
                    m_2,
                    eps_denom,
                    Difference(),
                )
            end
        end
    end

    return dephasing_time_min
end