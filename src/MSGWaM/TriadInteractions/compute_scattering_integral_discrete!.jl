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
    (; wavespectrum, col_int, partition) = spec_tend
    (; rhobar) = state.atmosphere

    was = @ivy view(wavespectrum, ii, jj, kk, :, :)
    sqrtrhobar = sqrt(rhobar[ii, jj, kk]) ##background density at the level (ii, jj, kk)


    # No update_interpolation_coef! here:
    # k is discrete and interpolation is only performed in m.

    dkp = kp[2] - kp[1]

    # Fourier-mode numbers represented by the spectral grid.
    nmin = round(Int, kp[1] / dkp)
    nmax = nmin + kpl - 1

    ntasks = state.namelists.triad.nthreads_triad

    @sync for tid in 1:ntasks
        inds = partition[tid]

        @spawn begin
            @ivy for idx in inds

                mi = (idx - 1) ÷ kpl + 1
                kpi = (idx - 1) % kpl + 1

                nk = was[kpi, mi]

                kr = kp[kpi]
                mr = m[mi]

                # Fourier-mode number corresponding to kr.
                nr = nmin + kpi - 1

                stk_sum = 0.0
                stk_diff = 0.0
                # ==================================================
                # Sum interactions
                #
                #       nr = n1 + n2
                # ==================================================

                if nr >= 2 * nmin
                    for n1 in nmin:(nr - nmin)

                        n2 = nr - n1

                        if n2 < nmin || n2 > nmax
                            continue
                        end

                        kp1i = n1 - nmin + 1
                        kp2i = n2 - nmin + 1

                        stk_sum += compute_st_k(
                            spec_tend, was, kp1i, kp2i, nk, kr, mr, triad_mode, Sum()
                        )
                    end
                end

                # ==================================================
                # Difference interactions
                #
                #       n1 = nr + n2
                # ==================================================

                if nr + nmin <= nmax
                    for n2 in nmin:(nmax - nr)

                        n1 = nr + n2

                        kp1i = n1 - nmin + 1
                        kp2i = n2 - nmin + 1

                        stk_diff += compute_st_k(
                            spec_tend, was, kp1i, kp2i, nk, kr, mr,
                            triad_mode, Difference()
                        )
                    end
                end

                # For the discrete spectrum N_n(m), the Δk factors
                # cancel. No horizontal quadrature weight is needed.

                col_int[ii, jj, kk, kpi, mi] =
                    4π * (stk_sum - stk_diff) / sqrtrhobar
            end
        end 
    end

    return nothing
end