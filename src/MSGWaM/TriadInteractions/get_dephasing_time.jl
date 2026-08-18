function get_dephasing_time end

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
    (; wavespectrum, col_int, action_ref) = spec_tend
    (; kp, m, delkp, delm, kpl) = spec_tend.spec_grid

    (; n2) = state.atmosphere
    (; u) = state.variables.predictands

    # No action-active spectral cell has a nonzero collision
    # integral in this physical grid cell.
    if isinf(tau_nl)
        return Inf
    end

    # Minimum significant spectral-cell action.
    action_floor = action_rel_tol * action_ref[]

    # Maximum nonlinear rate among action-active spectral cells.
    max_rate = 1.0 / tau_nl

    # Retain a parent spectral cell when
    #
    #     (|St| / N) / max_rate > increment_rel_tol.
    rate_cutoff = increment_rel_tol * max_rate

    n_local = sqrt(n2[ii, jj, kk])
    dudz = compute_dphidz_center(u, state, ii, jj, kk, identity)
    dndz = compute_dphidz_center(n2, state, ii, jj, kk, sqrt)

    dephasing_time_min = Inf
    eps_denom = 1.0e-14

    # For the continuous-k formulation only.
    if x_size > 1
        (; aa, la, qq, lq) = spec_tend.kin_box
    end

    # For the discrete-k formulation, determine the Fourier-mode
    # numbers represented by kp.
    if x_size == 1
        dkp = kp[2] - kp[1]
        nmin = round(Int, kp[1] / dkp)
        nmax = nmin + kpl - 1
    end

    for mi in eachindex(m), kpi in eachindex(kp)

        was = wavespectrum[ii, jj, kk, kpi, mi]
        st = col_int[ii, jj, kk, kpi, mi]

        #------------------------------------------------------
        # Contained wave action of the parent spectral cell.
        #
        # x_size > 1:
        #     A = N Δkp Δm
        #
        # x_size == 1:
        #     A = N_k Δm
        #------------------------------------------------------

        spectral_cell_measure = if x_size == 1
            abs(delm[mi])
        else
            abs(delkp[kpi] * delm[mi])
        end

        cell_action = was * spectral_cell_measure

        # Exclude spectral cells whose contained action is below
        # the initialized action floor.
        if cell_action <= action_floor
            continue
        end

        # was must be positive here because cell_action > action_floor.
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

        # ======================================================
        # DISCRETE HORIZONTAL WAVENUMBERS
        # ======================================================

        if x_size == 1

            n_parent = nmin + kpi - 1

            #--------------------------------------------------
            # Sum interactions:
            #
            #     n_parent = n1 + n2
            #--------------------------------------------------

            if n_parent >= 2 * nmin

                for n1 in nmin:(n_parent - nmin)

                    n2_mode = n_parent - n1

                    if n2_mode < nmin || n2_mode > nmax
                        continue
                    end

                    kp1i = n1 - nmin + 1
                    kp2i = n2_mode - nmin + 1

                    kp_1 = kp[kp1i]
                    kp_2 = kp[kp2i]

                    # Positive resonance branch.
                    m_1, m_2 = compute_m1m2(
                        kp_parent, kp_1, kp_2, m_parent, Sum(), Sum()
                    )

                    if check_resolved_spectral_mode(spec_tend, kp_1, m_1, Triad2D()) &&
                       check_resolved_spectral_mode(spec_tend, kp_2, m_2, Triad2D())

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

                    # Negative resonance branch.
                    m_1, m_2 = compute_m1m2(
                        kp_parent, kp_1, kp_2, m_parent, Sum(), Difference()
                    )

                    if check_resolved_spectral_mode(spec_tend, kp_1, m_1, Triad2D()) &&
                       check_resolved_spectral_mode(spec_tend, kp_2, m_2, Triad2D())

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
            end

            #--------------------------------------------------
            # Difference interactions:
            #
            #     n1 = n_parent + n2
            #--------------------------------------------------

            if n_parent + nmin <= nmax

                for n2_mode in nmin:(nmax - n_parent)

                    n1 = n_parent + n2_mode

                    kp1i = n1 - nmin + 1
                    kp2i = n2_mode - nmin + 1

                    kp_1 = kp[kp1i]
                    kp_2 = kp[kp2i]

                    # Positive resonance branch.
                    m_1, m_2 = compute_m1m2(
                        kp_parent,
                        kp_1,
                        kp_2,
                        m_parent,
                        Difference(),
                        Sum(),
                    )

                    if check_resolved_spectral_mode(spec_tend, kp_1, m_1, Triad2D()) &&
                       check_resolved_spectral_mode(spec_tend, kp_2, m_2, Triad2D())

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

                    # Negative resonance branch.
                    m_1, m_2 = compute_m1m2(
                        kp_parent,
                        kp_1,
                        kp_2,
                        m_parent,
                        Difference(),
                        Difference(),
                    )

                    if check_resolved_spectral_mode(spec_tend, kp_1, m_1, Triad2D()) &&
                       check_resolved_spectral_mode(spec_tend, kp_2, m_2, Triad2D())

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

        # ======================================================
        # CONTINUOUS HORIZONTAL WAVENUMBERS
        # ======================================================

        else

            aar = aa[kpi]
            qqr = qq[kpi]

            #--------------------------------------------------
            # Sum interaction manifold
            #--------------------------------------------------

            for ai in 1:la[kpi]

                #----------------------------------------------
                # Left branch
                #----------------------------------------------

                p_left = aar[ai] - kp_parent
                kp_1, kp_2 = compute_kp1kp2(kp_parent, p_left, Sum())

                m_1, m_2 = compute_m1m2(
                    kp_parent, kp_1, kp_2, m_parent, Sum(), Sum()
                )

                if check_resolved_spectral_mode(spec_tend, kp_1, m_1, Triad2D()) &&
                   check_resolved_spectral_mode(spec_tend, kp_2, m_2, Triad2D())

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

                m_1, m_2 = compute_m1m2(
                    kp_parent, kp_1, kp_2, m_parent, Sum(), Difference()
                )

                if check_resolved_spectral_mode(spec_tend, kp_1, m_1, Triad2D()) &&
                   check_resolved_spectral_mode(spec_tend, kp_2, m_2, Triad2D())

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

                #----------------------------------------------
                # Right branch
                #----------------------------------------------

                p_right = kp_parent - aar[ai]
                kp_1, kp_2 = compute_kp1kp2(kp_parent, p_right, Sum())

                m_1, m_2 = compute_m1m2(
                    kp_parent, kp_1, kp_2, m_parent, Sum(), Sum()
                )

                if check_resolved_spectral_mode(spec_tend, kp_1, m_1, Triad2D()) &&
                   check_resolved_spectral_mode(spec_tend, kp_2, m_2, Triad2D())

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

                m_1, m_2 = compute_m1m2(
                    kp_parent, kp_1, kp_2, m_parent, Sum(), Difference()
                )

                if check_resolved_spectral_mode(spec_tend, kp_1, m_1, Triad2D()) &&
                   check_resolved_spectral_mode(spec_tend, kp_2, m_2, Triad2D())

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

            # --------------------------------------------------
            # Difference interaction manifold
            # --------------------------------------------------

            for qi in 1:lq[kpi]

                q_val = qqr[qi]
                kp_1, kp_2 = compute_kp1kp2(
                    kp_parent, q_val, Difference()
                )

                m_1, m_2 = compute_m1m2(
                    kp_parent,
                    kp_1,
                    kp_2,
                    m_parent,
                    Difference(),
                    Sum(),
                )

                if check_resolved_spectral_mode(spec_tend, kp_1, m_1, Triad2D()) &&
                   check_resolved_spectral_mode(spec_tend, kp_2, m_2, Triad2D())

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

                m_1, m_2 = compute_m1m2(
                    kp_parent,
                    kp_1,
                    kp_2,
                    m_parent,
                    Difference(),
                    Difference(),
                )

                if check_resolved_spectral_mode(spec_tend, kp_1, m_1, Triad2D()) &&
                   check_resolved_spectral_mode(spec_tend, kp_2, m_2, Triad2D())

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
    end

    return dephasing_time_min
end