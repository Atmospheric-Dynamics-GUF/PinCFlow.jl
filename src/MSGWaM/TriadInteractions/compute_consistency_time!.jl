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
    (; kp, m) = spec_tend.spec_grid
    (; aa, la, qq, lq) = spec_tend.kin_box

    (; n2) = state.atmosphere
    (; u) = state.variables.predictands
    (; wavespectrum, col_int) = spec_tend
    mu_pl = 1.0
    eps_denom = 1.0e-14

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1

        n_local = sqrt(n2[i, j, k])
        dudz = compute_dphidz_center(u, state, i, j, k, identity)
        dndz = compute_dphidz_center(n2, state, i, j, k, sqrt)

        consistency_time_min = Inf

        for kpi_parent in eachindex(kp), mi_parent in eachindex(m)
            if wavespectrum[i, j, k, kpi_parent, mi_parent] > eps_denom && abs(col_int[i, j, k, kpi_parent, mi_parent]) > eps_denom
                kp_parent = kp[kpi_parent]
                m_parent = m[mi_parent]

                abs(m_parent) < eps_denom && continue

                aar = aa[kpi_parent]
                qqr = qq[kpi_parent]

                # -------------------------
                # Sum interaction manifold
                # -------------------------
                for ii in 1:la[kpi_parent]

                    # left branch
                    p_left = aar[ii] - kp_parent

                    kp_1_res, kp_2_res =
                        compute_kp1kp2(kp_parent, p_left, Sum())

                    kpi_1 = compute_nearest_index(kp, kp_1_res)
                    kpi_2 = compute_nearest_index(kp, kp_2_res)

                    m_1_res, m_2_res =
                        compute_m1m2(kp_parent, kp_1_res, kp_2_res, m_parent, Sum(), Sum())

                    mi_1 = compute_nearest_index(m, m_1_res)
                    mi_2 = compute_nearest_index(m, m_2_res)

                    consistency_time_min =
                        update_consistency_time(
                            consistency_time_min,
                            n_local,
                            dudz,
                            dndz,
                            kp,
                            m,
                            kpi_parent,
                            mi_parent,
                            kpi_1,
                            mi_1,
                            kpi_2,
                            mi_2,
                            eps_denom,
                            Sum(),
                        )

                    m_1_res, m_2_res =
                        compute_m1m2(kp_parent, kp_1_res, kp_2_res, m_parent, Sum(), Difference())

                    mi_1 = compute_nearest_index(m, m_1_res)
                    mi_2 = compute_nearest_index(m, m_2_res)

                    consistency_time_min =
                        update_consistency_time(
                            consistency_time_min,
                            n_local,
                            dudz,
                            dndz,
                            kp,
                            m,
                            kpi_parent,
                            mi_parent,
                            kpi_1,
                            mi_1,
                            kpi_2,
                            mi_2,
                            eps_denom,
                            Sum(),
                        )

                    # right branch
                    p_right = kp_parent - aar[ii]

                    kp_1_res, kp_2_res =
                        compute_kp1kp2(kp_parent, p_right, Sum())

                    kpi_1 = compute_nearest_index(kp, kp_1_res)
                    kpi_2 = compute_nearest_index(kp, kp_2_res)

                    m_1_res, m_2_res =
                        compute_m1m2(kp_parent, kp_1_res, kp_2_res, m_parent, Sum(), Sum())

                    mi_1 = compute_nearest_index(m, m_1_res)
                    mi_2 = compute_nearest_index(m, m_2_res)

                    consistency_time_min =
                        update_consistency_time(
                            consistency_time_min,
                            n_local,
                            dudz,
                            dndz,
                            kp,
                            m,
                            kpi_parent,
                            mi_parent,
                            kpi_1,
                            mi_1,
                            kpi_2,
                            mi_2,
                            eps_denom,
                            Sum(),
                        )

                    m_1_res, m_2_res =
                        compute_m1m2(kp_parent, kp_1_res, kp_2_res, m_parent, Sum(), Difference())

                    mi_1 = compute_nearest_index(m, m_1_res)
                    mi_2 = compute_nearest_index(m, m_2_res)

                    consistency_time_min =
                        update_consistency_time(
                            consistency_time_min,
                            n_local,
                            dudz,
                            dndz,
                            kp,
                            m,
                            kpi_parent,
                            mi_parent,
                            kpi_1,
                            mi_1,
                            kpi_2,
                            mi_2,
                            eps_denom,
                            Sum(),
                        )
                end

                # -------------------------------
                # Difference interaction manifold
                # -------------------------------
                for jj in 1:lq[kpi_parent]

                    q_val = qqr[jj]

                    kp_1_res, kp_2_res =
                        compute_kp1kp2(kp_parent, q_val, Difference())

                    kpi_1 = compute_nearest_index(kp, kp_1_res)
                    kpi_2 = compute_nearest_index(kp, kp_2_res)

                    m_1_res, m_2_res =
                        compute_m1m2(kp_parent, kp_1_res, kp_2_res, m_parent, Difference(), Sum())

                    mi_1 = compute_nearest_index(m, m_1_res)
                    mi_2 = compute_nearest_index(m, m_2_res)

                    consistency_time_min =
                        update_consistency_time(
                            consistency_time_min,
                            n_local,
                            dudz,
                            dndz,
                            kp,
                            m,
                            kpi_parent,
                            mi_parent,
                            kpi_1,
                            mi_1,
                            kpi_2,
                            mi_2,
                            eps_denom,
                            Difference(),
                        )

                    m_1_res, m_2_res =
                        compute_m1m2(kp_parent, kp_1_res, kp_2_res, m_parent, Difference(), Difference())

                    mi_1 = compute_nearest_index(m, m_1_res)
                    mi_2 = compute_nearest_index(m, m_2_res)

                    consistency_time_min =
                        update_consistency_time(
                            consistency_time_min,
                            n_local,
                            dudz,
                            dndz,
                            kp,
                            m,
                            kpi_parent,
                            mi_parent,
                            kpi_1,
                            mi_1,
                            kpi_2,
                            mi_2,
                            eps_denom,
                            Difference(),
                        )
                end
            end
        end

        consistency_time[i, j, k] = mu_pl * consistency_time_min
    end

    return nothing
end