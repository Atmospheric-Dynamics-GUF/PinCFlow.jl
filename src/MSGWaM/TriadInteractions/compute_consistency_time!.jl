function compute_consistency_time! end

function compute_consistency_time!(state::State)
    (; wkb_mode) = state.namelists.wkb
    (; triad_mode) = state.namelists.triad
    compute_consistency_time!(state, wkb_mode, triad_mode)
    return
end

function compute_consistency_time!(
    state::State,
    wkb_mode::Union{MultiColumn, SingleColumn},
    triad_mode::Triad2D;
)
    (; domain) = state
    (; i0, i1, j0, j1, k0, k1) = domain

    (; spec_tend) = state.wkb
    (; consistency_time) = spec_tend
    (; kp, m, delkp, delm) = spec_tend.spec_grid
    (; aa, la, qq, lq) = spec_tend.kin_box

    (; n2) = state.atmosphere
    (; u) = state.variables.predictands
    (; wavespectrum, col_int) = spec_tend
    (; action_abs_tol, st_abs_tol) = state.namelists.triad

    eps_denom = 1.0e-14

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1

        n_local = sqrt(n2[i, j, k])
        dudz = compute_dphidz_center(u, state, i, j, k, identity)
        dndz = compute_dphidz_center(n2, state, i, j, k, sqrt)

        max_was = maximum(wavespectrum[i, j, k, :, :])

        if max_was <= 1.0E-15
            #return if there is non significant wad in the physical grid cell
            continue
        end

        consistency_time_min = Inf

        for mi in eachindex(m), kpi in eachindex(kp)
            dkps = delkp[kpi]
            dms = delm[mi]
            wave_action = wavespectrum[i, j, k, kpi, mi] * dkps * dms
            if wave_action > action_abs_tol && abs(col_int[i, j, k, kpi, mi]) > st_abs_tol
                kp_parent = kp[kpi]
                m_parent = m[mi]

                abs(m_parent) < eps_denom && continue

                aar = aa[kpi]
                qqr = qq[kpi]

                # -------------------------
                # Sum interaction manifold
                # -------------------------
                for ii in 1:la[kpi]

                    # left branch
                    p_left = aar[ii] - kp_parent

                    kp_1, kp_2 =
                        compute_kp1kp2(kp_parent, p_left, Sum())

                    m_1, m_2 =
                        compute_m1m2(kp_parent, kp_1, kp_2, m_parent, Sum(), Sum())

                    consistency_time_min =
                        update_consistency_time(
                            consistency_time_min,
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
                            Sum()
                        )

                    m_1, m_2 =
                        compute_m1m2(kp_parent, kp_1, kp_2, m_parent, Sum(), Difference())

                    consistency_time_min =
                        update_consistency_time(
                            consistency_time_min,
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
                            Sum()
                        )
                    # right branch
                    p_right = kp_parent - aar[ii]

                    kp_1, kp_2 =
                        compute_kp1kp2(kp_parent, p_right, Sum())

                    m_1, m_2 =
                        compute_m1m2(kp_parent, kp_1, kp_2, m_parent, Sum(), Sum())


                    consistency_time_min =
                        update_consistency_time(
                            consistency_time_min,
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
                            Sum()
                        )

                    m_1, m_2 =
                        compute_m1m2(kp_parent, kp_1, kp_2, m_parent, Sum(), Difference())

                    consistency_time_min =
                        update_consistency_time(
                            consistency_time_min,
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
                            Sum()
                        )
                end

                # -------------------------------
                # Difference interaction manifold
                # -------------------------------
                for jj in 1:lq[kpi]

                    q_val = qqr[jj]

                    kp_1, kp_2 =
                        compute_kp1kp2(kp_parent, q_val, Difference())

                    m_1, m_2 =
                        compute_m1m2(kp_parent, kp_1, kp_2, m_parent, Difference(), Sum())


                    consistency_time_min =
                        update_consistency_time(
                            consistency_time_min,
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
                            Difference()
                        )

                    m_1, m_2 =
                        compute_m1m2(kp_parent, kp_1, kp_2, m_parent, Difference(), Difference())


                    consistency_time_min =
                        update_consistency_time(
                            consistency_time_min,
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
                            Difference()
                        )
                end
            end
        end
        consistency_time[i, j, k] =  consistency_time_min
    end

    return nothing
end