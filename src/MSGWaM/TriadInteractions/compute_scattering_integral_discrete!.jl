function compute_scattering_integral_discrete! end

function compute_scattering_integral_discrete!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    triad_mode::Triad2D,
)
    (; spec_tend) = state
    (; kp, m, kpl) = spec_tend.spec_grid
    (; wavespectrum, col_int, diag_dephasing_time, partition) = spec_tend
    (; n2, rhobar) = state.atmosphere
    (; u) = state.variables.predictands
    (; compute_dephasing_time, nthreads_triad) = state.namelists.triad

    was = @ivy view(wavespectrum, ii, jj, kk, :, :)

    rhobar_local = rhobar[ii, jj, kk]

    if compute_dephasing_time
        n_local = sqrt(n2[ii, jj, kk])
        dudz = compute_dphidz_center(u, state, ii, jj, kk, identity)
        dndz = compute_dphidz_center(n2, state, ii, jj, kk, sqrt)
    else
        n_local = 0.0
        dudz = 0.0
        dndz = 0.0
    end

    # No update_interpolation_coef! here:
    # k is discrete and interpolation is only performed in m.

    dkp = kp[2] - kp[1]

    nmin = round(Int, kp[1] / dkp)
    nmax = nmin + kpl - 1

    @sync for tid in 1:nthreads_triad
        inds = partition[tid]

        @spawn begin
            @ivy for idx in inds
                mi = (idx - 1) ÷ kpl + 1
                kpi = (idx - 1) % kpl + 1

                nk = was[kpi, mi]

                kr = kp[kpi]
                mr = m[mi]

                nr = nmin + kpi - 1

                stk_sum = 0.0
                stk_diff = 0.0
                dephasing_time_k = Inf

                # ==========================================================
                # Sum interactions
                #
                #     nr = n1 + n2
                # ==========================================================

                if nr >= 2 * nmin
                    for n1 in nmin:(nr - nmin)
                        n2_mode = nr - n1

                        if n2_mode < nmin || n2_mode > nmax
                            continue
                        end

                        kp1i = n1 - nmin + 1
                        kp2i = n2_mode - nmin + 1

                        st_value, tau_dep = compute_st_k(spec_tend, was, kp1i, kp2i, nk, kr, mr, n_local, dudz, dndz, compute_dephasing_time, triad_mode, Sum())

                        stk_sum += st_value
                        dephasing_time_k = min(dephasing_time_k, tau_dep)
                    end
                end

                # ==========================================================
                # Difference interactions
                #
                #     n1 = nr + n2
                # ==========================================================

                if nr + nmin <= nmax
                    for n2_mode in nmin:(nmax - nr)
                        n1 = nr + n2_mode

                        kp1i = n1 - nmin + 1
                        kp2i = n2_mode - nmin + 1

                        st_value, tau_dep = compute_st_k(spec_tend, was, kp1i, kp2i, nk, kr, mr, n_local, dudz, dndz, compute_dephasing_time, triad_mode, Difference())

                        stk_diff += st_value
                        dephasing_time_k = min(dephasing_time_k, tau_dep)
                    end
                end

                # For the discrete spectrum N_n(m), the Δk factors cancel.
                # No horizontal quadrature weight is required.

                col_int[ii, jj, kk, kpi, mi] =
                    4π * (stk_sum - stk_diff) / rhobar_local

                diag_dephasing_time[kpi, mi] = dephasing_time_k
            end
        end
    end

    return nothing
end